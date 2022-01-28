/*!
Pass a custom match and traceback function to the WFA algorithm

Constraints:
 - max qlen and tlen is max value of i32 because of
   * hi and lo fields in WaveFrontSet
   * further bound by computing a_k in [wf_align]

 */
use num;
use std::cmp;

mod utils;
pub mod types;



pub fn wf_traceback(
    all_wavefronts: &types::WaveFronts,
    score: usize,
    config: &types::Config
) -> String {
    if config.verbosity > 1 {
        eprintln!("\t[wfa::wf_backtrace]");
    }

    let x: i32 = config.penalties.mismatch;
    let o: i32 = config.penalties.gap_open;
    let e: i32 = config.penalties.gap_extend;

    let mut cigar = String::new();

    // start with central diagonal
    let mut k = all_wavefronts.a_k as i32;

    // start at the furthest offset on the m-wavefront i.e. the end of the alignment
    let m_wf =  all_wavefronts.get_m_wavefront(score).unwrap();
    let wave_length = m_wf.len();
    let hi = m_wf.hi;
    let k_index = utils::compute_k_index(wave_length, k, hi);

    // offset
    let m_s_k = m_wf.offsets[k_index];

    let mut v: usize = utils::compute_v(m_s_k, k);
    let mut h: usize = utils::compute_h(m_s_k, k);

    if config.verbosity >  5 {
        eprintln!("\t\t({}, {})", v, h);
    }

    let mut backtrace_op = types::BacktraceOperation::MatchMismatch;

    let mut s = score as i32;
    let mut offset = m_s_k as i32;

    while v > 0 && h > 0 && s > 0 {
        // compute scores
        let gap_open_score: i32 = s - o - e;
        let gap_extend_score: i32 = s - e;
        let mismatch_score: i32 = s - x;

        if config.verbosity >  5 {
            eprintln!("\t\tscore: {} \n\
                       \t\tOperation: {:?} \n\
                       \t\t{{\n\
                       \t\t\tg_o: {} \n\
                       \t\t\tg_e: {} \n\
                       \t\t\tx: {} \n\
                       \t\t}}\
                       ",
                      s, backtrace_op,
                      gap_open_score, gap_extend_score, mismatch_score);
        }

        let del_ext: Option<i32> = if backtrace_op == types::BacktraceOperation::Insertion {
            None
        } else {
            utils::backtrace_utils::backtrace_deletion_extend_offset(all_wavefronts, gap_extend_score, k)
        };

        let del_open: Option<i32> = if backtrace_op == types::BacktraceOperation::Insertion {
            None
        } else {
            utils::backtrace_utils::backtrace_deletion_open_offset(all_wavefronts, gap_open_score, k)
        };

        let ins_ext: Option<i32> = if backtrace_op == types::BacktraceOperation::Deletion {
            None
        } else {
            utils::backtrace_utils::backtrace_insertion_extend_offset(all_wavefronts, gap_extend_score, k)
        };

        let ins_open: Option<i32> = if backtrace_op == types::BacktraceOperation::Deletion {
            None
        } else {
            utils::backtrace_utils::backtrace_insertion_open_offset(all_wavefronts, gap_open_score, k)
        };

        let misms: Option<i32> = if backtrace_op != types::BacktraceOperation::MatchMismatch {
            None
        } else {
            utils::backtrace_utils::backtrace_mismatch_offset(all_wavefronts, mismatch_score, k)
        };

        // Compute maximum offset
        let max_all: Option<i32> = vec![del_ext, del_open, ins_ext, ins_open, misms]
            .into_iter()
            .max()
            .unwrap();

        if config.verbosity > 5 {
            let res = vec![del_ext, del_open, ins_ext, ins_open, misms];
            eprintln!("\t\tdel_ext, del_open, ins_ext, ins_open, misms\n\
                       \t\tops {:?} \n\
                       \t\toffset {} \n\
                       \t\tmax_all {:?} \n\
                       \t\tbacktrace_op {:?}",
                      res, offset, max_all, backtrace_op);
        }

        // Traceback Matches
        if max_all.is_some() &&
            backtrace_op == types::BacktraceOperation::MatchMismatch &&
            offset >= max_all.unwrap() {
            let num_matches = (offset - max_all.unwrap()) as u32;
            utils::backtrace_utils::backtrace_matches_check(
                &mut offset,
                &mut cigar,
                num_matches,
                k
            );
            offset = max_all.unwrap();
        }

        if max_all == del_ext {
            // Extend a deletion
            cigar.push('D');
            // Update state
            s = gap_extend_score;
            k += 1;
            backtrace_op = types::BacktraceOperation::Deletion;
        } else if max_all == del_open {
            // Open a deletion
            cigar.push('D');
            // Update state
            s = gap_open_score;
            k += 1;
            backtrace_op = types::BacktraceOperation::MatchMismatch;
        } else if max_all == ins_ext {
            // Extend an insertion
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

        if config.verbosity > 5 {
            eprintln!("\t\t({}, {}) s {}", v, h, s);
        }
    }

    // reached the end of one or both of the sequences
    if s == 0 {
        // backtrace matches check
        let num_matches = offset as u32;
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

    let reversed_cigar = cigar.chars().rev().collect::<String>();
    reversed_cigar
}

pub fn wf_extend<F>(
    m_wavefront: &mut types::WaveFront,
    match_lambda: &F,
    config: &types::Config
)
where
    F: Fn(usize, usize) -> bool,
{
    if config.verbosity > 1 {
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

            if config.verbosity > 254 {
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

    let hi_or_min = |maybe_wf: Option<&types::WaveFront>| -> i32 {
        match maybe_wf {
            Some(wf) => wf.hi,
            _ => i32::MIN
        }
    };

    let lo_or_max = |maybe_wf: Option<&types::WaveFront>| -> i32 {
        match maybe_wf {
            Some(wf) => wf.lo,
            _ => i32::MAX
        }
    };

    let hi: i32 = *vec![
        hi_or_min(wavefronts.get_m_wavefront(utils::to_usize_or_zero(s_x))),
        hi_or_min(wavefronts.get_m_wavefront(utils::to_usize_or_zero(s_o_e))),
        hi_or_min(wavefronts.get_i_wavefront(utils::to_usize_or_zero(s_e))),
        hi_or_min(wavefronts.get_d_wavefront(utils::to_usize_or_zero(s_e))),
    ].iter().max().unwrap() as i32 + 1;

    let lo: i32 = *vec![
        lo_or_max(wavefronts.get_m_wavefront(utils::to_usize_or_zero(s_x))),
        lo_or_max(wavefronts.get_m_wavefront(utils::to_usize_or_zero(s_o_e))),
        lo_or_max(wavefronts.get_i_wavefront(utils::to_usize_or_zero(s_e))),
        lo_or_max(wavefronts.get_d_wavefront(utils::to_usize_or_zero(s_e))),
    ].iter().min().unwrap() - 1;

    (hi, lo)
}

pub fn wf_next(
    wavefronts: &mut types::WaveFronts,
    score: usize,
    config: &types::Config
) {
    if config.verbosity > 1 {
        eprintln!("\t[wflambda::wf_next]");
    }

    if config.verbosity > 5 {
        eprintln!("\t\tscore {}", score);
    }

    let s: i32 = score as i32;
    let x: i32 = config.penalties.mismatch;
    let o: i32 = config.penalties.gap_open;
    let e: i32 = config.penalties.gap_extend;

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
        if config.verbosity > 5 {
            eprintln!("\t\tskipping score {}", score);
        }
            return;
    }

    if config.verbosity > 5 {
        eprintln!("\t\ts {} s - o - e {} s - e {} s - x {}",
                  s, unsigned_s_o_e, unsigned_s_e, unsigned_s_x);
        eprint!("\t\tk\tI\tD\tM");
        eprintln!();
    }

    // compute the highest/rightmost and lowest/leftmost diagonal for a
    // wavefront with the given score will reach
    let (hi, lo) = compute_wf_next_limits(&wavefronts, score as isize);

    // Vec::<types::WfType>::new();
    let wavefronts_to_allocate = vec![types::WfType::I, types::WfType::D, types::WfType::M];

    // Allocate the next wave front
    let allocation_result = wavefronts.allocate_wavefronts(score as u32, lo, hi, &wavefronts_to_allocate);
    match allocation_result {
        Ok(_) => {
            if config.verbosity > 254 {
                eprintln!("\t\tallocated wfs for score: {} \n\
                           \t {:?}"
                          , score, wavefronts.wavefront_set[score])
            }
        },
        Err(msg) => { panic!("{}", msg) }
    };

    let wave_length: usize = utils::compute_wave_length(lo, hi);

    for k in lo..=hi {

        // eprintln!("\t\twavelength {}", wave_length);

        if k-1 < lo || k+1 > hi {
            continue;
        }

        let k_index_sub_one: usize = utils::compute_k_index(wave_length, k-1, hi);
        let k_index: usize = utils::compute_k_index(wave_length, k, hi);
        let k_index_add_one: usize = utils::compute_k_index(wave_length, k+1, hi);

        if  config.verbosity > 254 {
            eprintln!("\t\t\tk: {}\tk-1 {}\tk-index {}\tk+1 {}",
                      k, k_index_sub_one, k_index, k_index_add_one);
        }

        let i_s_k: Option<i32> = vec![
            wavefronts
                .get_m_wavefront(unsigned_s_o_e)
                .and_then(|m_wf| m_wf.get_offset(k-1))
                .cloned(),
            wavefronts
                .get_i_wavefront(unsigned_s_e)
                .and_then(|i_wf| i_wf.get_offset(k-1))
                .cloned(),
        ].into_iter().max().unwrap().map(|x| x + 1);

        let d_s_k: Option<i32> = vec![
            wavefronts
                .get_m_wavefront(unsigned_s_o_e)
                .and_then(|m_wf| m_wf.get_offset(k+1))
                .cloned(),
            wavefronts
                .get_d_wavefront(unsigned_s_e)
                .and_then(|d_wf| d_wf.get_offset(k+1))
                .cloned(),
        ].into_iter().max().unwrap();

        let m_s_k: Option<i32> = vec![
            wavefronts
                .get_m_wavefront(unsigned_s_x)
                .and_then(|m_wf| m_wf.get_offset(k))
                .cloned()
                .map(|x| x+1),
            i_s_k,
            d_s_k,
        ].into_iter().max().unwrap();

        // offsets
        if config.verbosity > 5 {
            eprint!("\t\t{}\t{:?}\t{:?}\t{:?}", k, i_s_k, d_s_k, m_s_k);
            eprintln!();
        }


        // set the values
        let wf_set: &mut Option<types::WaveFrontSet> = &mut wavefronts.wavefront_set[score];
        // eprintln!("\t\t 1 {:?}", wf_set);
        // eprintln!();
        let wf_set: &mut types::WaveFrontSet = wf_set.as_mut().unwrap();
        // eprintln!("\t\t 2 {:?}", wf_set);

        match i_s_k {
            Some(o) => {
                match  wf_set.i.as_mut() {
                    Some(i) => {i.offsets[k_index] = o}
                    _ => {}
                }
            },
            _ => { wf_set.i = None }
        };

        match d_s_k {
            Some(o) => {
                match  wf_set.d.as_mut() {
                    Some(d) => {d.offsets[k_index] = o}
                    _ => {}
                }
            },

            _ => { wf_set.d = None }
        };
        match m_s_k {
            Some(o) => {
                match  wf_set.m.as_mut() {
                    Some(m) => {m.offsets[k_index] = o}
                    _ => {}
                }
            },
            _ => { wf_set.m = None }
        };
    }

    if config.verbosity > 4 {
        eprintln!();
    }
}


pub fn wf_align<F>(
    tlen: u32,
    qlen: u32,
    match_lambda: &F,
    config: &types::Config
) -> (usize, String)
where
    F: Fn(usize, usize) -> bool
{
    if config.verbosity > 1 {
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
        i: None,
        d: None,
        m: Some(types::WaveFront::new(hi, lo)),
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
    if config.verbosity > 0 {
        eprintln!("Config\n\
                   \ttlen: {}\n\
                   \tqlen: {}\n\
                   \ta_k: {}\n\
                   \ta_offset: {}",
                  qlen, tlen, a_k, a_offset);
    }

    if config.verbosity > 2 {
        eprintln!("Loop");
    }

    loop {
        // Extend the current wavefront
        if all_wavefronts.get_m_wavefront(score).is_some() {
            let m_wf_mut: &mut types::WaveFront = &mut all_wavefronts
                .wavefront_set[score]
                .as_mut()
                .unwrap()
                .m
                .as_mut()
                .unwrap();
            wf_extend(m_wf_mut, match_lambda, &config);
        }

        // Check whether we have reached the final point
        // Get the m-wavefront with the current score
        if utils::end_reached(all_wavefronts.get_m_wavefront(score), a_k, a_offset) {
            if config.verbosity > 1 {
                let m_s: &types::WaveFront = all_wavefronts.get_m_wavefront(score).unwrap();
                let i_s: &types::WaveFront = all_wavefronts.get_m_wavefront(score).unwrap();
                let d_s: &types::WaveFront = all_wavefronts.get_m_wavefront(score).unwrap();

                eprintln!("\n--------------------------------------\n");
                eprintln!("\twf\tscore\tk\toffset");

                for (name, wf) in vec![("M", m_s), ("I", i_s), ("D", d_s)].iter() {
                    for k in wf.lo..=wf.hi {
                        let k_index = utils::compute_k_index(wf.len(), k as i32, wf.hi);
                        let offset = wf.offsets[k_index];
                        eprintln!("\t{}\t{}\t{}\t{}",
                                  name, score, k, offset);
                    }
                    eprintln!();
                }

                eprintln!("\n--------------------------------------\n");
            }
            break
        }

        score += 1;

        // compute the next wavefront
        wf_next(&mut all_wavefronts, score, config);
    }

    // let cigar = String::new();
    let cigar = wf_traceback(&all_wavefronts, score, config);

    let each_wf = vec![ types::WfType::M ];
    utils::debug_utils::visualize_all(&all_wavefronts, a_offset*2, &each_wf);

    (score, cigar)
}


#[cfg(test)]
mod tests {
    mod test_config {
        pub static CONFIG: crate::types::Config = crate::types::Config {
            adapt: false,
            verbosity: 0,
            penalties: crate::types:: Penalties {
                mismatch: 4,
                matches: 0,
                gap_open: 6,
                gap_extend: 2,
            },
        };
    }

    mod core_functions {
        #[test]
        fn test_wf_next() {
            assert!(false);
        }
    }

    mod align {

        use super::{super::*, *};

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

            let (score, cigar) = wf_align(tlen as u32, qlen as u32, &match_lambda, &test_config::CONFIG);

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

                let (score, cigar) = wf_align(tlen  as u32, qlen  as u32, &match_lambda, &test_config::CONFIG);
                eprintln!("Result:\n\tScore: {} Cigar {}", score, cigar);
                crate::utils::backtrace_utils::print_aln(&cigar[..], t, q);
            }

            {
                // different sequences
                let text  = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\
                             TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                             XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
                let query = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\
                             TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT\
                             XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";

                let tlen = text.len();
                let qlen = query.len();

                let t: &[u8] = text.as_bytes();
                let q: &[u8] = query.as_bytes();

                let match_lambda = |v: usize, h: usize| {
                    v < tlen && h < qlen && t[v] == q[h]
                };

                // let (score, cigar) = wf_align(tlen as u32, qlen as u32, &match_lambda, &test_config::CONFIG);
                // eprintln!("Result:\n\tScore: {} Cigar {}", score, cigar);
                // crate::utils::backtrace_utils::print_aln(&cigar[..], t, q);
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

                // let (score, cigar) = wf_align(tlen as u32, qlen as u32, &match_lambda, &test_config::CONFIG);
                // eprintln!("Result:\n\tScore: {} Cigar {}", score, cigar);
                // crate::utils::backtrace_utils::print_aln(&cigar[..], t, q);

            }
        }
    }
}
