/*!
Pass a custom match and traceback function to the WFA algorithm

Constraints:
 - max qlen and tlen is max value of i32 because of
   * hi and lo fields in WaveFrontSet
   * further bound by computing a_k in [wf_align]

 */
mod utils;
pub mod types;
mod config;

use num;
use std::cmp;

const NULL_OFFSETT: i32 = 0;

/// Change isizes to i32
pub fn wf_traceback(
		all_wavefronts: &types::WaveFronts,
		score: usize
) -> String {
		let NULL_OFFSET = -10;
		if config::VERBOSITY > 1 {
        eprintln!("\t[wfa::wf_backtrace]");
    }

		let x: isize = 4;
    let o: isize = 6;
    let e: isize = 2;

		let mut cigar = String::new();

		// start with central diagonal
		let mut k = all_wavefronts.a_k as i32;

		let m_wf = &all_wavefronts.wavefront_set[score].m;
		let len = m_wf.vals.len();
		let hi = m_wf.hi;
		let k_index = utils::compute_k_index(len, k, hi);

		// offset
		let m_s_k = m_wf.vals[k_index];

		let mut v: usize = utils::compute_v(m_s_k, k);
		let mut h: usize = utils::compute_h(m_s_k, k);

		let mut backtrace_op = types::BacktraceOperation::MatchMismatch;

		let mut s = score as isize;
    let mut offset = m_s_k as isize;

		while v > 0 && h > 0 && s > 0 {
				// compute scores
        let gap_open_score = s - o - e;
        let gap_extend_score = s - e;
        let mismatch_score = s - x;



				let del_ext: isize = if backtrace_op == types::BacktraceOperation::Insertion {
						NULL_OFFSET
				} else {
						utils::backtrace_utils::backtrace_deletion_extend_offset(all_wavefronts, s, k).unwrap_or(NULL_OFFSET)
				};

				let del_open: isize = if backtrace_op == types::BacktraceOperation::Insertion {
            NULL_OFFSET
        } else {
            utils::backtrace_utils::backtrace_deletion_open_offset(all_wavefronts, gap_open_score, k).unwrap_or(NULL_OFFSET)
        };

        let ins_ext: isize = if backtrace_op == types::BacktraceOperation::Deletion {
            NULL_OFFSET
        } else {
            utils::backtrace_utils::backtrace_insertion_extend_offset(all_wavefronts, gap_extend_score, k).unwrap_or(NULL_OFFSET)
        };

        let ins_open: isize = if backtrace_op == types::BacktraceOperation::Deletion {
            NULL_OFFSET
        } else {
            utils::backtrace_utils::backtrace_insertion_open_offset(all_wavefronts, gap_open_score, k).unwrap_or(NULL_OFFSET)
        };

        let misms: isize = if backtrace_op != types::BacktraceOperation::MatchMismatch {
            NULL_OFFSET
        } else {
            utils::backtrace_utils::backtrace_mismatch_offset(all_wavefronts, mismatch_score, k).unwrap_or(NULL_OFFSET)
        };

				// Compute maximum offset
        let max_all: isize = *vec![del_ext, del_open, ins_ext, ins_open, misms]
						.iter()
						.max()
						.unwrap();

				// Traceback Matches
        if backtrace_op == types::BacktraceOperation::MatchMismatch && offset >= max_all {
            let num_matches = (offset - max_all) as usize;
            utils::backtrace_utils::backtrace_matches_check(
                &mut offset,
                &mut cigar,
                num_matches,
                k
            );
            offset = max_all;
        }

				if max_all == del_ext {
            // Add Deletion
            cigar.push('D');
            // Update state
            s = gap_extend_score;
            k += 1;
            backtrace_op = types::BacktraceOperation::Deletion;
        } else if max_all == del_open {
            // Add Deletion
            cigar.push('D');
            // Update state
            s = gap_open_score;
            k += 1;
            backtrace_op = types::BacktraceOperation::MatchMismatch;
        } else if max_all == ins_ext {
            // Add Insertion
            cigar.push('I');
            // Update state
            s = gap_extend_score;
            k -= 1;
            offset -= 1;
            backtrace_op = types::BacktraceOperation::Insertion;
        } else if max_all == ins_open {
            // Add Insertion
            cigar.push('I');
            // Update state
            s = gap_open_score;
            k -= 1;
            offset -= 1;
            backtrace_op = types::BacktraceOperation::MatchMismatch;
        } else if max_all == misms {
            // Add Mismatch
            cigar.push('X');

            // Update state
            s = mismatch_score;
            offset -= 1;
        } else {
						panic!("Backtrace error: No link found during backtrace");
        }

				v = crate::utils::compute_v(offset as i32, k);
				h = crate::utils::compute_h(offset as i32, k);
		};

		// reached the end of one or both of the sequences
		if score == 0 {
				// backtrace matches check
        let num_matches = offset as usize;
				utils::backtrace_utils::backtrace_matches_check(
            &mut offset,
            &mut cigar,
            num_matches,
            k
        );
		} else {
				// add indels
				while v > 0 {
            cigar.push('D');
            v -= 1;
        }

        while h > 0 {
            cigar.push('I');
            h -= 1;
        }
		}

		cigar
}

pub fn wf_extend<F>(
		m_wavefront: &mut types::WaveFront,
		match_lambda: &F,
)
where
		F: Fn(usize, usize) -> bool,
{
		if config::VERBOSITY > 1 {
				eprintln!("\t[wflambda::wf_extend]");
		}

		// eprintln!("\t\tlo {} hi {}",  m_wavefront.lo, m_wavefront.hi);

		for k in m_wavefront.lo..=m_wavefront.hi {

				let k_index: usize = utils::compute_k_index(m_wavefront.vals.len(), k, m_wavefront.hi);

				// assuming tlen > qlen
				let m_s_k: i32 = m_wavefront.vals[k_index];

				let mut v: usize = utils::compute_v(m_s_k, k);
				let mut h: usize = utils::compute_h(m_s_k, k);

				// eprintln!("\t\t\tk {}\toffset {}\t({}, {})", k, m_s_k, v, h);
				// vt[v][h] = m_wavefront.vals[k_index] as i32;

				while match_lambda(v, h) {

						if config::VERBOSITY > 254 {
								eprintln!("\t[wflambda::wf_extend]\n\
													 \t\t({} {})",
													v, h);
						}

						// increment our offset on the m_wavefront
						m_wavefront.vals[k_index] += 1;

						v+=1;
						h+=1;
				}
		}
}

fn compute_wf_next_limits(
		wavefronts: &types::WaveFronts,
		score: isize
) -> (i32, i32) {
		let s: isize = score;
		let x: isize = 4;
		let o: isize = 6;
		let e: isize = 2;

		let s_x: isize = s - x;
		let s_o_e: isize = s - o - e;
		let s_e: isize = s - e;

		let hi: i32 = *vec![
				wavefronts.get(utils::to_usize_or_zero(s_x)).m.hi,
				wavefronts.get(utils::to_usize_or_zero(s_o_e)).m.hi,
				wavefronts.get(utils::to_usize_or_zero(s_e)).i.hi,
				wavefronts.get(utils::to_usize_or_zero(s_e)).d.hi,
		].iter().max().unwrap() as i32 + 1;

		let lo: i32 = *vec![
				wavefronts.get(utils::to_usize_or_zero(s_x)).m.lo,
				wavefronts.get(utils::to_usize_or_zero(s_o_e)).m.lo,
				wavefronts.get(utils::to_usize_or_zero(s_e)).i.lo,
				wavefronts.get(utils::to_usize_or_zero(s_e)).d.lo,
		].iter().min().unwrap() - 1;

		(hi, lo)
}

pub fn wf_next(
		wavefronts: &mut types::WaveFronts,
		score: usize
) {
		if config::VERBOSITY > 1 {
				eprintln!("\t[wflambda::wf_next]");
		}

		// compute the highest/rightmost and lowest/leftmost diagonal for a
		// wavefront with the given score will reach
		let (hi, lo) = compute_wf_next_limits(&wavefronts, score as isize);

		// Allocate the next wave front
		wavefronts.allocate_wavefronts(score as u32, lo, hi).unwrap();

		let s: i32 = score as i32;
		let x: i32 = 4;
		let o: i32 = 6;
		let e: i32 = 2;

		let signed_s_x: i32 = s - x;
		let signed_s_o_e: i32 = s - o - e;
		let signed_s_e: i32 = s - e;

		if config::VERBOSITY > 4 {
				eprintln!("\t\tscore {}\n\
									 \t\tsigned_s_x {} signed_s_o_e {} signed_s_e {}\n\
									 \t\tlo: {} hi: {}\n",
									score, signed_s_x, signed_s_o_e, signed_s_e, lo, hi);
		}

		let s_x = signed_s_x as usize;
		let s_o_e = signed_s_o_e as usize;
		let s_e = signed_s_e as usize;

		if config::VERBOSITY > 254 {
				eprintln!("\t\tk-1\tk\tk+1");
		}

		if config::VERBOSITY > 4 {
				eprintln!("\t\tscore {}\n\
									 \t\ts_x {} s_o_e {} s_e {}\n\
									 \t\tlo: {} hi: {}\n",
									score, s_x, s_o_e, s_e, lo, hi);
		}

		if config::VERBOSITY > 4 {
				eprint!("\t\tk\tI\tD\tM");
				eprintln!();
		}

		let wave_length: usize = utils::compute_wave_length(lo, hi);

		for k in lo..=hi {

				// eprintln!("\t\twavelength {}", wave_length);

				if k-1 < lo || k+1 > hi {
						continue;
				}

				let k_index_sub_one: usize = utils::compute_k_index(wave_length, k-1, hi);
				let k_index: usize = utils::compute_k_index(wave_length, k, hi);
				let k_index_add_one: usize = utils::compute_k_index(wave_length, k+1, hi);

				if  config::VERBOSITY > 254 {
						eprintln!("\t\t{}\t{}\t{}", k_index_add_one, k_index, k_index_add_one);
				}

				let mut i_s_k: i32 = *vec![
						match wavefronts.option_get(s_o_e) {
								Some(wf) => {
										wf.m.vals.get(k_index_sub_one).unwrap_or(&signed_s_o_e)
								},
								_ => &signed_s_o_e
						},
						match wavefronts.option_get(s_e) {
								Some(wf) => wf.i.vals.get(k_index_sub_one).unwrap_or(&signed_s_e),
								_ => &signed_s_e
						},
				].iter().max().unwrap() + 1;

				let d_s_k: i32 = **vec![
						match wavefronts.option_get(s_o_e) {
								Some(wf) => wf.m.vals.get(k_index_add_one).unwrap_or(&signed_s_o_e),
								_ => &NULL_OFFSETT
						},
						match wavefronts.option_get(s_e) {
								Some(wf) => wf.d.vals.get(k_index_add_one).unwrap_or(&signed_s_e),
								_ => &signed_s_e
						},
				].iter().max().unwrap();

				let m_s_k: i32 = *vec![
						match wavefronts.option_get(s_x) {
								Some(wf) => *wf.m.vals.get(k_index).unwrap_or(&signed_s_x) + 1,
								_ => NULL_OFFSETT + 1
						},
						i_s_k,
						d_s_k,
				].iter().max().unwrap();

				// offsets
				if config::VERBOSITY > 4 {
						eprint!("\t\t{}\t{}\t{}\t{}", k, i_s_k, d_s_k, m_s_k);
						eprintln!();
				}

				// set the values
				wavefronts.wavefront_set[score].i.vals[k_index] = i_s_k;
				wavefronts.wavefront_set[score].d.vals[k_index] = d_s_k;
				wavefronts.wavefront_set[score].m.vals[k_index] = m_s_k;
		}

		if config::VERBOSITY > 4 {
				eprintln!();
		}
}

pub fn wf_align<F>(
		tlen: u32,
		qlen: u32,
		match_lambda: &F
) -> (usize, String)
where
		F: Fn(usize, usize) -> bool
{
		if config::VERBOSITY > 1 {
				eprintln!("[wflambda::wf_align]");
		}

		// compute the central diagonal, a_k.
		let a_k: usize =  num::abs_sub(tlen as isize, qlen as isize) as usize;

		// the furthest offset we expect the central diagonal to reach
		// subtract 1 because of the zero index
		let a_offset: u32 = cmp::max(tlen, qlen);

		// eprintln!("\t a_k {} a_offset {}", a_k, a_offset);

		let hi: i32 = a_k as i32;
		let lo: i32 = a_k as i32;

		// Initial conditions
		let wf_set = types::WaveFrontSet {
				i: types::WaveFront::new(hi, lo),
				d: types::WaveFront::new(hi, lo),
				m: types::WaveFront::new(hi, lo),
		};

		let mut all_wavefronts = types::WaveFronts {
				wavefront_set: vec![wf_set],
				min_k: -(cmp::min(tlen, qlen) as i32),
				max_k: cmp::max(tlen, qlen) as i32,
				a_k,
		};

		// score
		let mut score: usize = 0;

		// set score at start to 0
		// unnecessary
		// score ... diagonal
		all_wavefronts.wavefront_set[0].m.vals[0] = 0;

		// Print config
		if config::VERBOSITY > 0 {
				eprintln!("Config\n\
									 \ttlen: {}\n\
									 \tqlen: {}\n\
									 \ta_k: {}\n\
									 \ta_offset: {}",
									qlen, tlen, a_k, a_offset);
		}

		if config::VERBOSITY > 2 {
				eprintln!("Loop");
		}

		loop {
				// Extend the current wavefront
				wf_extend(&mut all_wavefronts.wavefront_set[score].m, match_lambda);

				// Check whether we have reached the final point
				// Get the m-wavefront with the current score
				let m_wavefront: &types::WaveFront = &all_wavefronts.wavefront_set[score].m;
				let k_index: usize = utils::compute_k_index(m_wavefront.vals.len(), a_k as i32, m_wavefront.hi);

				if config::VERBOSITY > u8::MAX - 1 {
				eprintln!("\n\
									 \t[wflambda::wf_align]\n\
									 \t\tlo:{} hi: {} len: {} k_index: {}",
									m_wavefront.lo, m_wavefront.hi,
									m_wavefront.vals.len(), k_index);
				}

				// current offset on the m wavefront at the current score
				let m_s_k = m_wavefront.vals[k_index];

				if m_s_k >= 0 && (m_s_k as u32) >= a_offset {
						eprintln!("Done");
						eprintln!("----------------------");
						eprintln!("m_s_k {}, score {}", m_s_k, score);

						eprintln!("----------------------");

						break
				}

				if config::VERBOSITY > 2 {
						eprintln!("\n-------------------------------\n\
											 \t[wflambda::wf_align]\n\
											 \t\tscore: {}\n\
											 \t\toffset: {}\n",
											score, m_s_k);
				}

				score += 1;

				// compute the next wavefront
				wf_next(&mut all_wavefronts, score);
		}

		let cigar = String::new();
		let cigar = wf_traceback(&all_wavefronts, score);

		let each_wf = vec![ types::WfType::M ];
		utils::debug_utils::visualize_all(&all_wavefronts, a_offset, &each_wf);

		(score, cigar)
}


#[cfg(test)]
mod tests {

		mod core_functions {
				#[test]
				fn test_wf_next() {
						assert!(false);
				}
		}

		mod align {

				use super::super::*;

				#[test]
				fn align_same_sequence() {
						// different sequences
						let text  = "GAGAAT";
						let query = "GAGAAT";

						let tlen = text.len();
						let qlen = query.len();

						let t: &[u8] = text.as_bytes();
						let q: &[u8] = query.as_bytes();

						let match_lambda = |v: usize, h: usize| {
								v < tlen && h < qlen && t[v] == q[h]
						};

						let (score, cigar) = wf_align(tlen as u32, qlen as u32, &match_lambda);

						assert_eq!(score, 0);

						eprintln!("\nScore: {}", score);
						eprintln!();
						utils::backtrace_utils::print_aln(&cigar[..], t, q);
				}

				#[test]
				fn align_different_sequence_same_len() {
						{
								// different sequences
								let text =  "GAGATA";
								let query = "GACACA";

								let tlen = text.len();
								let qlen = query.len();

								let t: &[u8] = text.as_bytes();
								let q: &[u8] = query.as_bytes();

								let match_lambda = |v: usize, h: usize| {
										v < tlen && h < qlen && t[v] == q[h]
								};

								eprintln!("--------------------");

								let (score, cigar) = wf_align(tlen  as u32, qlen  as u32, &match_lambda);

								eprintln!("--------------------");

								eprintln!("Result:\n\tScore: {}", score);

								eprintln!("--------------------");

								utils::backtrace_utils::print_aln(&cigar[..], t, q);
						}

						{
								// different sequences
								let text =  "GAGATA";
								let query = "GACACA";

								let tlen = text.len();
								let qlen = query.len();

								let t: &[u8] = text.as_bytes();
								let q: &[u8] = query.as_bytes();

								let match_lambda = |v: usize, h: usize| {
										v < tlen && h < qlen && t[v] == q[h]
								};

								// let (score, cigar) = wf_align(tlen, qlen, &match_lambda);
								// assert_eq!(score, 8);
								// utils::backtrace_utils::print_aln(&cigar[..], t, q);
						}

						{
								// different sequences
								let text  = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
								let query = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";

								let tlen = text.len();
								let qlen = query.len();

								let t: &[u8] = text.as_bytes();
								let q: &[u8] = query.as_bytes();

								let match_lambda = |v: usize, h: usize| {
										v < tlen && h < qlen && t[v] == q[h]
								};

								// let score = wf_align(tlen as u32, qlen as u32, &match_lambda);
						}

						{
								// different sequences
								let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
								let text  = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";

								let tlen = text.len();
								let qlen = query.len();

								let t: &[u8] = text.as_bytes();
								let q: &[u8] = query.as_bytes();

								let match_lambda = |v: usize, h: usize| {
										v < tlen && h < qlen && t[v] == q[h]
								};

								// let (score, cigar) = wf_align(tlen as u32, qlen as u32, &match_lambda);
						}
				}
		}
}
