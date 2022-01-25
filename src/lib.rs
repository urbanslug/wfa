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

				let k_index: usize = utils::compute_k_index(m_wavefront.len(), k, m_wavefront.hi);

				// assuming tlen > qlen
				let m_s_k: i32 = m_wavefront.offsets[k_index];

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
						m_wavefront.offsets[k_index] += 1;

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
				wavefronts.get_m_wavefront(utils::to_usize_or_zero(s_x)).unwrap().hi,
				wavefronts.get_m_wavefront(utils::to_usize_or_zero(s_o_e)).unwrap().hi,
				wavefronts.get_i_wavefront(utils::to_usize_or_zero(s_e)).unwrap().hi,
				wavefronts.get_d_wavefront(utils::to_usize_or_zero(s_e)).unwrap().hi,
		].iter().max().unwrap() as i32 + 1;

		let lo: i32 = *vec![
				wavefronts.get_m_wavefront(utils::to_usize_or_zero(s_x)).unwrap().lo,
				wavefronts.get_m_wavefront(utils::to_usize_or_zero(s_o_e)).unwrap().lo,
				wavefronts.get_i_wavefront(utils::to_usize_or_zero(s_e)).unwrap().lo,
				wavefronts.get_d_wavefront(utils::to_usize_or_zero(s_e)).unwrap().lo,
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

		if config::VERBOSITY > 2 {
				eprintln!("\t\tscore {}", score);
		}

		let s: i32 = score as i32;
		let x: i32 = 4;
		let o: i32 = 6;
		let e: i32 = 2;

		let signed_s_x: i32 = s - x;
    let signed_s_o_e: i32 = s - o - e;
    let signed_s_e: i32 = s - e;


    let unsigned_s_x = signed_s_x as usize;
    let unsigned_s_o_e = signed_s_o_e as usize;
    let unsigned_s_e = signed_s_e as usize;

    if (signed_s_x < 0 || wavefronts.get_m_wavefront(unsigned_s_x).is_none()) &&
        (signed_s_o_e < 0 || wavefronts.get_m_wavefront(unsigned_s_o_e).is_none()) &&
        (signed_s_e < 0 || wavefronts.get_i_wavefront(unsigned_s_e).is_none()) &&
        (signed_s_e < 0 ||  wavefronts.get_d_wavefront(unsigned_s_e).is_none())
		{

						eprintln!("skipping score {}", score);
						return;
		}

		/*
		if config::VERBOSITY > 4 {
				eprintln!("\t\tscore {}\n\
									 \t\tsigned_s_x {} signed_s_o_e {} signed_s_e {}\n\
									 \t\tlo: {} hi: {}\n",
									score, signed_s_x, signed_s_o_e, signed_s_e, lo, hi);
		}
		 */

		if config::VERBOSITY > 254 {
				eprintln!("\t\tk-1\tk\tk+1");
		}

		/*
		if config::VERBOSITY > 4 {
				eprintln!("\t\tscore {}\n\
									 \t\ts_x {} s_o_e {} s_e {}\n\
									 \t\tlo: {} hi: {}\n",
									score, s_x, s_o_e, s_e, lo, hi);
		}
		 */

		if config::VERBOSITY > 4 {
				eprint!("\t\tk\tI\tD\tM");
				eprintln!();
		}

		// compute the highest/rightmost and lowest/leftmost diagonal for a
		// wavefront with the given score will reach
		let (hi, lo) = compute_wf_next_limits(&wavefronts, score as isize);

		// Allocate the next wave front
		wavefronts.allocate_wavefronts(score as u32, lo, hi).unwrap();

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

				let i_s_k: i32 = *vec![
						match wavefronts.get_m_wavefront(unsigned_s_o_e) {
								Some(m_wf) => m_wf.offsets.get(k_index_sub_one).unwrap_or(&signed_s_o_e),
								_ => &signed_s_o_e
						},
						match wavefronts.get_i_wavefront(unsigned_s_e) {
								Some(i_wf) => i_wf.offsets.get(k_index_sub_one).unwrap_or(&signed_s_e),
								_ => &signed_s_e
						},
				].iter().max().unwrap() + 1;

				let d_s_k: i32 = **vec![
						match wavefronts.get_m_wavefront(unsigned_s_o_e) {
								Some(m_wf) => m_wf.offsets.get(k_index_add_one).unwrap_or(&signed_s_o_e),
								_ => &signed_s_o_e
						},
						match wavefronts.get_d_wavefront(unsigned_s_e) {
								Some(d_wf) => d_wf.offsets.get(k_index_add_one).unwrap_or(&signed_s_e),
								_ => &signed_s_e
						},

				].iter().max().unwrap();

				let m_s_k: i32 = *vec![
						match wavefronts.get_m_wavefront(unsigned_s_x) {
								Some(m_wf) => *m_wf.offsets.get(k_index).unwrap_or(&signed_s_x) + 1,
								_ => &signed_s_x + 1
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
				let wf_set: &mut Option<types::WaveFrontSet> = &mut wavefronts.wavefront_set[score];
				let wf_set: &mut types::WaveFrontSet = wf_set.as_mut().unwrap();
				wf_set.i.offsets[k_index] = i_s_k;
				wf_set.d.offsets[k_index] = d_s_k;
				wf_set.m.offsets[k_index] = m_s_k;
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
				wavefront_set: vec![Some(wf_set)],
				min_k: -(cmp::min(tlen, qlen) as i32),
				max_k: cmp::max(tlen, qlen) as i32,
				a_k,
		};

		// score
		let mut score: usize = 0;

		// set score at start to 0
		// unnecessary
		// score ... diagonal
		// all_wavefronts.wavefront_set[0].m.vals[0] = 0;
		assert_eq!(all_wavefronts.get_m_wavefront(score).unwrap().offsets[a_k], 0);

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
				if all_wavefronts.get_m_wavefront(score).is_some() {
						let m_wf_mut: &mut types::WaveFront = &mut all_wavefronts
								.wavefront_set[score]
								.as_mut()
								.unwrap()
								.m;
						wf_extend(m_wf_mut, match_lambda);
				}

				// Check whether we have reached the final point
				// Get the m-wavefront with the current score
				if utils::end_reached(all_wavefronts.get_m_wavefront(score), a_k, a_offset) {
						break
				}

				score += 1;

				// compute the next wavefront
				wf_next(&mut all_wavefronts, score);
		}

		let cigar = String::new();
		// let cigar = wf_traceback(&all_wavefronts, score);

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
