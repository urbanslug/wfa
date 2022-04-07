/*!
Base level WFA
 */
use num;
use std::cmp;

use crate::types;
use crate::utils::{self, backtrace as backtrace_utils};
use crate::core;

fn wf_traceback(
    all_wavefronts: &types::WaveFronts,
    score: usize,
    config: &types::Config,
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
    let m_wf = all_wavefronts.get_m_wavefront(score as i32).unwrap();
    let wave_length = m_wf.len();
    let hi = m_wf.hi;
    let k_index = utils::compute_k_index(wave_length, k, hi);

    // offset
    let m_s_k = m_wf.offsets[k_index];

    let mut v = utils::compute_v(m_s_k, k);
    let mut h = utils::compute_h(m_s_k, k);

    if config.verbosity > 5 {
        eprintln!("\t\t({}, {})", v, h);
    }

    let mut backtrace_op = types::BacktraceOperation::MatchMismatch;

    let mut s = score as i32;
    let mut offset = m_s_k;

    while v > 0 && h > 0 && s > 0 {
        // compute scores
        let gap_open_score: i32 = s - o - e;
        let gap_extend_score: i32 = s - e;
        let mismatch_score: i32 = s - x;

        if config.verbosity > 4 {
            eprintln!(
                "\t\tscore: {} \n\
                       \t\tOperation: {:?} \n\
                       \t\t{{\n\
                       \t\t\tg_o: {} \n\
                       \t\t\tg_e: {} \n\
                       \t\t\tx: {} \n\
                       \t\t}}\
                       ",
                s, backtrace_op, gap_open_score, gap_extend_score, mismatch_score
            );
        }

        let del_ext: Option<i32> = if backtrace_op == types::BacktraceOperation::Insertion {
            None
        } else {
            backtrace_utils::backtrace_deletion_extend_offset(
                all_wavefronts,
                gap_extend_score,
                k,
            )
        };

        let del_open: Option<i32> = if backtrace_op == types::BacktraceOperation::Insertion {
            None
        } else {
            backtrace_utils::backtrace_deletion_open_offset(
                all_wavefronts,
                gap_open_score,
                k,
            )
        };

        let ins_ext: Option<i32> = if backtrace_op == types::BacktraceOperation::Deletion {
            None
        } else {
            backtrace_utils::backtrace_insertion_extend_offset(
                all_wavefronts,
                gap_extend_score,
                k,
            )
        };

        let ins_open: Option<i32> = if backtrace_op == types::BacktraceOperation::Deletion {
            None
        } else {
            backtrace_utils::backtrace_insertion_open_offset(
                all_wavefronts,
                gap_open_score,
                k,
            )
        };

        let misms: Option<i32> = if backtrace_op != types::BacktraceOperation::MatchMismatch {
            None
        } else {
            backtrace_utils::backtrace_mismatch_offset(all_wavefronts, mismatch_score, k)
        };

        // Compute maximum offset
        let max_all: Option<i32> = vec![del_ext, del_open, ins_ext, ins_open, misms]
            .into_iter()
            .max()
            .unwrap();

        if config.verbosity > 4 {
            let res = vec![del_ext, del_open, ins_ext, ins_open, misms];
            eprintln!(
                "\t\tdel_ext, del_open, ins_ext, ins_open, misms\n\
                       \t\tops {:?} \n\
                       \t\toffset {} \n\
                       \t\tmax_all {:?} \n\
                       \t\tbacktrace_op {:?}",
                res, offset, max_all, backtrace_op
            );
        }

        // Traceback Matches
        if max_all.is_some()
            && backtrace_op == types::BacktraceOperation::MatchMismatch
            && offset >= max_all.unwrap()
        {
            let num_matches = (offset - max_all.unwrap()) as u32;
            backtrace_utils::backtrace_matches_check(
                &mut offset,
                &mut cigar,
                num_matches);
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

        v = crate::utils::compute_v(offset, k);
        h = crate::utils::compute_h(offset, k);

        if config.verbosity > 5 {
            eprintln!("\t\t({}, {}) s {}", v, h, s);
        }
    }

    // reached the end of one or both of the sequences
    if s == 0 {
        // backtrace matches check
        let num_matches = offset as u32;
        backtrace_utils::backtrace_matches_check(&mut offset, &mut cigar, num_matches);
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

// TODO: remove match_positions debug_vec
fn wf_extend(
    m_wavefront: &mut types::WaveFront,
    config: &types::Config,
    match_positions: &mut Vec<(i32, i32, usize)>,
    _score: usize,
    text: &[u8],
    query: &[u8],
) {
    if config.verbosity > 1 {
        eprintln!("\t[wfa::wf_extend]");
    }

    let tlen = text.len();
    let qlen = query.len();


    // eprintln!("\t\tlo {} hi {}",  m_wavefront.lo, m_wavefront.hi);
   //  eprintln!("\t\tscore {}", score);

    for k in m_wavefront.lo..=m_wavefront.hi {

        let m_s_k: i32 = *m_wavefront.get_offset(k).expect("[wfa::wf_extend] fail unwrap k={k}");
        // let k_index: usize = utils::new_compute_k_index(k, m_wavefront.lo, m_wavefront.hi);

        // assuming tlen > qlen
        // let m_s_k: i32 = m_wavefront.offsets[k_index];

        let mut v = utils::compute_v(m_s_k, k);
        let mut h = utils::compute_h(m_s_k, k);

        eprintln!("\t\t\tk {}\toffset {}\t({}, {})", k, m_s_k, v, h);
        // vt[v][h] = m_wavefront.vals[k_index] as i32;

        let is_match = |v: i32, h: i32| -> bool {
            if v < 0 || h < 0 || h  as usize >= tlen || v as usize >= qlen {
                return false;
            }
            let v = v as usize;
            let h = h as usize;

            eprintln!("\t\t\t {} {}", text[h] as char, query[v] as char);
            text[h] == query[v]
        };

        while is_match(v, h) {
            if config.verbosity > 254 {
                eprintln!(
                    "\t[wflambda::wf_extend]\n\
                           \t\t({} {})",
                    v, h
                );
            }
            // increment our offset on the m_wavefront
            *m_wavefront.get_offset_mut(k).unwrap() += 1;
            // m_wavefront.offsets[k_index] += 1;

            // Debug stuff
            {
                let offset = *m_wavefront
                    .get_offset(k)
                    .expect("DEBUG [wfa::wf_extend] fail unwrap k={k}");
                match_positions.push((v, h, offset as usize));

                eprintln!("\t\tk {}\toffset {}", k, offset);
            }

            v += 1;
            h += 1;
        }
    }
}

pub fn wf_align(
    text: &[u8],
    query: &[u8],
    config: &types::Config,
) -> (usize, String) {
    if config.verbosity > 1 {
        eprintln!("[wfa::wf_align]");
    }

    let tlen = text.len() as u32;
    let qlen = query.len() as u32;

    // compute the central diagonal, a_k.
    let a_k: i32 = (tlen as isize - qlen as isize) as i32;

    // the furthest offset we expect the central diagonal to reach
    // subtract 1 because of the zero index
    let a_offset: u32 = tlen;

    let hi: i32 = 0;
    let lo: i32 = 0;

    eprintln!("\t a_k {a_k} a_offset {a_offset} hi {hi} lo {lo}");

    // Initial conditions
    let wf_set = types::WaveFrontSet {
        i: None,
        d: None,
        m: Some(types::WaveFront::new(hi, lo)),
    };

    let mut all_wavefronts = types::WaveFronts {
        wavefront_set: vec![Some(wf_set)],
        min_k: -(qlen as isize),
        max_k: tlen as isize,
        a_k,
    };

    // score
    let mut score: usize = 0;

    // set score at start to 0
    // unnecessary
    // score ... diagonal
    // all_wavefronts.wavefront_set[0].m.vals[0] = 0;
    // Make sure start score is 0
    if *all_wavefronts
        .get_m_wavefront(score as i32)
        .expect(&format!("[wfa::wf_align] no m-wavefront at score {}", score))
        .get_offset(0 as i32)
        .expect(&format!("[wfa::wf_align] no offset on a_k = {a_k} m-wavefront at score {score}")) != 0
    {
            panic!("[wfa::wf_align] start score should be zero");
    }

    let mut match_posititons: Vec<(i32, i32, usize)> = Vec::new();

    // Print config
    if config.verbosity > 0 {
        eprintln!(
            "Config\n\
                   \ttlen: {}\n\
                   \tqlen: {}\n\
                   \ta_k: {}\n\
                   \ta_offset: {}",
            qlen, tlen, a_k, a_offset
        );
    }

    loop {
        // Extend the current wavefront
        if all_wavefronts.get_m_wavefront(score as i32).is_some() {
            let m_wf_mut: &mut types::WaveFront = &mut all_wavefronts.wavefront_set[score]
                .as_mut()
                .unwrap()
                .m
                .as_mut()
                .unwrap();

            wf_extend(
                m_wf_mut,
                &config,
                &mut match_posititons,
                score,
                text,
                query,
            );
        }

        // Check whether we have reached the final point
        // Get the m-wavefront with the current score
        if utils::end_reached(all_wavefronts.get_m_wavefront(score as i32), a_k, a_offset) {
            break;
        }

        score += 1;

        // compute the next wavefront
        core::wf_next(&mut all_wavefronts, score, config);
    }

    //let cigar = String::new();
    let cigar = wf_traceback(&all_wavefronts, score, config);

    (score, cigar)
}


#[cfg(test)]
mod tests {
    mod same_sequence {
        use crate::tests_prelude::*;

        #[test]
        fn test_short() {
            // different sequences
            let text  = "GAGAAT";
            let query = "GAGAAT";

            let t: &[u8] = text.as_bytes();
            let q: &[u8] = query.as_bytes();

            let (score, cigar) = wfa_align(t, q, &TEST_CONFIG);

            self::assert_eq!(score, 0);
            dbg!(score, &cigar);
        }
    }

    mod different_sequence_same_len {
        use crate::tests_prelude::*;

        #[test]
        fn test_short() {
            // different sequences
            let text = "GAGATA";
            let query = "GACACA";

            let t: &[u8] = text.as_bytes();
            let q: &[u8] = query.as_bytes();

            let (score, cigar) = wfa_align(t, q, &TEST_CONFIG);
            dbg!(score, &cigar);
            crate::utils::print_aln(cigar.as_bytes(), t, q);
        }

        #[test]
        fn test_long() {
            // different sequences
            let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
            let text = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";

            let t: &[u8] = text.as_bytes();
            let q: &[u8] = query.as_bytes();

            let (score, cigar) = wfa_align(t, q, &TEST_CONFIG);
            dbg!(score, &cigar);
            crate::utils::print_aln(cigar.as_bytes(), t, q);
        }

        #[test]
        fn test_longest() {
            // different sequences
            let text  = "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAA\
                         ATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGC\
                         GTTTGGAGAAATAAAATAGT";
            let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAA\
                         ATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCG\
                         CGTTGGAGAAATACAATAGT";

            let t: &[u8] = text.as_bytes();
            let q: &[u8] = query.as_bytes();

            let (score, cigar) = wfa_align(t, q, &TEST_CONFIG);
            dbg!(score, &cigar);
            crate::utils::print_aln(cigar.as_bytes(), t, q);
        }
    }

    mod different_sequence_different_len {
        use crate::tests_prelude::*;

        #[test]
        fn test_short_shorter_text() {
            // different sequences
            let text  =  "TCTGA";
            let query = "ATCTGA";

            // let text  =  "ACACA";
            // let query = "GACACA";

            let t: &[u8] = text.as_bytes();
            let q: &[u8] = query.as_bytes();

            let (score, cigar) = wfa_align(t, q, &TEST_CONFIG);
            dbg!(score, &cigar);

            self::assert_eq!(cigar, "DMMMMM");
            crate::utils::print_aln(cigar.as_bytes(), t, q);
        }

        #[test]
        fn test_short_shorter_query() {
            // different sequences
            let text  = "ATCTGA";
            let query =  "TCTGA";

            let t: &[u8] = text.as_bytes();
            let q: &[u8] = query.as_bytes();

            let (score, cigar) = wfa_align(t, q, &TEST_CONFIG);
            dbg!(score, &cigar);

            self::assert_eq!(cigar, "IMMMMM");
            crate::utils::print_aln(cigar.as_bytes(), t, q);

            // different sequences
            let text  = "ATCTGA";
            let query = "ATCGA";

            let t: &[u8] = text.as_bytes();
            let q: &[u8] = query.as_bytes();

            let (score, cigar) = wfa_align(t, q, &TEST_CONFIG);
            dbg!(score, &cigar);

            self::assert_eq!(cigar, "MMMIMM");
            crate::utils::print_aln(cigar.as_bytes(), t, q);
        }


        #[test]
        fn test_long_shorter_text() {
            // different sequences
            let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
            let text =  "TCTATACTGCGCGTTGAGAAATAAAAAG";

            let t: &[u8] = text.as_bytes();
            let q: &[u8] = query.as_bytes();

            let (score, cigar) = wfa_align(t, q, &TEST_CONFIG);

            self::assert_eq!(score, 40);
            crate::utils::print_aln(cigar.as_bytes(), t, q);
        }

        #[test]
        fn test_long_shorter_query() {
            // different sequences
            let query = "TCTTTATCGCGGTTGGAGAAAACAAAGT";
            let text  = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";

            let t: &[u8] = text.as_bytes();
            let q: &[u8] = query.as_bytes();

            let (score, cigar) = wfa_align(t, q, &TEST_CONFIG);
            self::assert_eq!(score, 40);
            crate::utils::print_aln(cigar.as_bytes(), t, q);
        }
    }
}
