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
    let m_wf =  all_wavefronts.get_m_wavefront(score as i32).unwrap();
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
    config: &types::Config,
    match_positions: &mut Vec<(usize, usize, usize)>
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

            let offset = m_wavefront.offsets[k_index];
            match_positions.push((v, h, offset as usize));

            v+=1;
            h+=1;
        }
    }
}

fn compute_wf_next_limits(
    wavefronts: &types::WaveFronts,
    score: usize,
    config: &types::Config
) -> (Option<i32>, Option<i32>) {
    enum WfLimit {
        Hi,
        Lo
    }

    let s: i32 = score as i32;
    let x: i32 = config.penalties.mismatch;
    let o: i32 = config.penalties.gap_open;
    let e: i32 = config.penalties.gap_extend;

    let s_x = s - x;
    let s_o_e = s - o - e;
    let s_e = s - e;

    let foo = |maybe_wf: Option<&types::WaveFront>, limit: WfLimit| -> Option<i32> {
        match maybe_wf {
            Some(wf) => match limit {
                WfLimit::Hi => Some(wf.hi),
                WfLimit::Lo => Some(wf.lo)
            },
            _ => None
        }
    };

    // what if a hi is empty?
    let hi: Option<i32> = vec![
        foo(wavefronts.get_m_wavefront(s_x), WfLimit::Hi),
        foo(wavefronts.get_m_wavefront(s_o_e), WfLimit::Hi),
        foo(wavefronts.get_i_wavefront(s_e), WfLimit::Hi),
        foo(wavefronts.get_d_wavefront(s_e), WfLimit::Hi),
    ].iter().max().unwrap().map(|x| x+1 );

    let maybe_los: Vec<Option<i32>> = vec![
        foo(wavefronts.get_m_wavefront(s_x), WfLimit::Lo),
        foo(wavefronts.get_m_wavefront(s_o_e), WfLimit::Lo),
        foo(wavefronts.get_i_wavefront(s_e), WfLimit::Lo),
        foo(wavefronts.get_d_wavefront(s_e), WfLimit::Lo),
    ].iter().filter(|x| x.is_some()).cloned().collect();

    if maybe_los.is_empty() {
        return (hi, None);
    }

    let lo: Option<i32> = maybe_los.iter().min().unwrap().map(|x| x-1);


    (hi, lo)
}

#[derive(Debug)]
struct AWFSet<'a>{
    // In
    in_m_sub: Option<&'a types::WaveFront>,
    in_m_gap: Option<&'a types::WaveFront>,
    in_i_ext: Option<&'a types::WaveFront>,
    in_d_ext: Option<&'a types::WaveFront>,

    // out
    out_m: Option<&'a types::WaveFront>,
    out_i: Option<&'a types::WaveFront>,
    out_d: Option<&'a types::WaveFront>,
}

fn fetch_wf<'a>(
    score: usize,
    wavefronts: &'a types::WaveFronts,
    config: &'a types::Config
) -> AWFSet<'a>
{
    let s: i32 = score as i32;
    let x: i32 = config.penalties.mismatch;
    let o: i32 = config.penalties.gap_open;
    let e: i32 = config.penalties.gap_extend;

    let s_x: i32 = s - x;
    let s_o_e: i32 = s - o - e;
    let s_e: i32 = s - e;

    let maybe_in_m_sub: Option<&types::WaveFront> = wavefronts.get_m_wavefront(s_x);
    let maybe_in_m_gap: Option<&types::WaveFront> = wavefronts.get_m_wavefront(s_o_e);
    let maybe_in_i_ext: Option<&types::WaveFront> = wavefronts.get_i_wavefront(s_e);
    let maybe_in_d_ext: Option<&types::WaveFront> = wavefronts.get_d_wavefront(s_e);

    AWFSet {
        in_m_sub: maybe_in_m_sub,
        in_m_gap: maybe_in_m_gap,
        in_i_ext: maybe_in_i_ext,
        in_d_ext: maybe_in_d_ext,

        out_m: None,
        out_i: None,
        out_d: None,
    }
}

fn alloc_wf(awf_set: &AWFSet) -> Vec<types::WfType> {
    let mut wavefronts_to_allocate = vec![types::WfType::M];

    // Allocate I-Wavefront
    if awf_set.in_m_gap.is_some() || awf_set.in_i_ext.is_some(){
        wavefronts_to_allocate.push(types::WfType::I);
    }

    // Allocate D-Wavefront
    if awf_set.in_m_gap.is_some() || awf_set.in_d_ext.is_some(){
        wavefronts_to_allocate.push(types::WfType::D);
    }

    wavefronts_to_allocate
}

pub fn wf_next(
    wavefronts: &mut types::WaveFronts,
    score: usize,
    config: &types::Config
) {
    if config.verbosity > 1 {
        eprintln!("\t[wflambda::wf_next]");
    }

    if config.verbosity > 4 {
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


    let cloned_wf = wavefronts.clone();
    let awf_set = fetch_wf(score, &cloned_wf, config);

    if awf_set.in_m_sub.is_none() &&
        awf_set.in_m_gap.is_none() &&
        awf_set.in_d_ext.is_none() &&
        awf_set.in_d_ext.is_none()
    {
        if config.verbosity > 4 {
            eprintln!("\t\tskipping score {}", score);
        }
        return;
    }

    if config.verbosity > 5 {
        eprintln!("\t\ts {} s - o - e {} s - e {} s - x {}", s, unsigned_s_o_e, unsigned_s_e, unsigned_s_x);
        eprint!("\t\tk\tI\tD\tM");
        eprintln!();
    }

    //----
    // Compute limits
    //----
    // compute the highest/rightmost and lowest/leftmost diagonal for a
    // wavefront with the given score will reach
    let (hi, lo): (Option<i32>, Option<i32>) = compute_wf_next_limits(&wavefronts, score, config);

    if hi.is_none() || lo.is_none() {
        if config.verbosity > 4 {
            eprintln!("\t\tskipping (invalid diagonals) hi {:?} lo {:?}", lo, hi);
        }
        return;
    }

    let (hi, lo) = (hi.unwrap(), lo.unwrap());

    eprintln!("\t\tscore: {} lo {} hi {}", score, lo, hi);

    // Allocate the next wave front
    let wavefronts_to_allocate = alloc_wf(&awf_set);
    // alloc_wf(&wavefronts, &awf_set, score, lo, hi);
    let allocation_result = wavefronts.allocate_wavefronts(score as u32, lo, hi, &wavefronts_to_allocate);
    match allocation_result {
        Ok(_) => {
            if config.verbosity > 254 {
                eprintln!("\t\tallocated wfs for score: {} \n\t {:?}", score, wavefronts.wavefront_set[score])
            }
        },
        Err(msg) => { panic!("{}", msg) }
    };

    /*
    // awf_set.out_i = wavefronts.get_i_wavefront(score as i32);
    // awf_set.out_d = wavefronts.get_d_wavefront(score as i32);
    // awf_set.out_m = wavefronts.get_m_wavefront(score as i32);
     */

    // let wave_length: usize = utils::compute_wave_length(lo, hi);

    let compute_k_index = |wavefront_length: usize, hi: i32, k: i32| -> usize {
        wavefront_length - ((hi - k) as usize) - 1
    };

    let compute_maybe_k_index = |wavefront_length: usize, hi: i32, k: i32| -> Option<usize> {
        let x = hi - k;
        if x < 0 || wavefront_length < (x as usize) || wavefront_length - (x as usize) < 1 {
            return None;
        }

        Some(wavefront_length - ((hi - k) as usize) - 1)
    };

    let assign_offsets_m = |wavefronts: &mut types::WaveFronts| {
        let wf_set: &mut Option<types::WaveFrontSet> = &mut wavefronts.wavefront_set[score];
        let wf_set: &mut types::WaveFrontSet = wf_set.as_mut().unwrap();
        let m_wf: &mut types::WaveFront = &mut wf_set.m.as_mut().unwrap();

        let in_m_wf: &types::WaveFront = awf_set.in_m_sub.unwrap();

        for k in lo..=hi {

            let k_index = compute_k_index(m_wf.len(), m_wf.hi, k);
            let maybe_in_k_index = compute_maybe_k_index(in_m_wf.len(), in_m_wf.hi, k);
            match maybe_in_k_index {
                Some(in_k_index) => {
                    m_wf.offsets[k_index] = in_m_wf.offsets[in_k_index] + 1;
                },
                None =>  {
                    m_wf.offsets[k_index] = -10;
                }
            };

            eprintln!("{} {}", k, m_wf.offsets[k_index]);
            eprintln!();
        }
    };

    let assign_offsets_im = |wavefronts: &mut types::WaveFronts| {
        let wf_set: &mut Option<types::WaveFrontSet> = &mut wavefronts.wavefront_set[score];
        let wf_set: &mut types::WaveFrontSet = wf_set.as_mut().unwrap();

        let out_m_wf: &mut types::WaveFront = &mut wf_set.m.as_mut().unwrap();
        let out_i_wf: &mut types::WaveFront = &mut wf_set.i.as_mut().unwrap();

        let in_m_sub_wf: &types::WaveFront = awf_set.in_m_sub.unwrap();
        let in_m_gap_wf: &types::WaveFront = awf_set.in_m_gap.unwrap();
        let in_i_ext_wf: &types::WaveFront = awf_set.in_i_ext.unwrap();

        for k in lo..=hi {
            // Update I
            let out_k_index: usize = compute_k_index(out_m_wf.len(), out_m_wf.hi, k);
            let maybe_in_k_index_gap: Option<usize> = compute_maybe_k_index(in_m_gap_wf.len(), in_m_gap_wf.hi, k-1);
            let maybe_in_k_index_ext: Option<usize> = compute_maybe_k_index(in_i_ext_wf.len(), in_i_ext_wf.hi, k-1);
            let x: Option<i32> = match maybe_in_k_index_gap {
                Some(in_k_index_gap) => Some(in_m_gap_wf.offsets[in_k_index_gap]),
                None => None,
            };

            let y: Option<i32> = match maybe_in_k_index_ext {
                Some(in_k_index_ext) => Some(in_i_ext_wf.offsets[in_k_index_ext]),
                None => None,
            };

            let maybe_ins: Option<i32> = *vec![x, y].iter().max().unwrap();
            match maybe_ins {
                Some(v) => out_i_wf.offsets[out_k_index] = v+1,
                _ => out_i_wf.offsets[out_k_index] = -10,
            }


            // Update M
            let out_k_index = compute_k_index(out_m_wf.len(), out_m_wf.hi, k);
            let maybe_in_k_index = compute_maybe_k_index(in_m_sub_wf.len(), in_m_sub_wf.hi, k);

            match maybe_in_k_index {
                Some(in_k_index) => {
                    out_m_wf.offsets[out_k_index] = in_m_sub_wf.offsets[in_k_index] + 1;
                },
                None =>  {
                    out_m_wf.offsets[out_k_index] = -10;
                }
            };
        }
    };

    let assign_offsets_dm = |wavefronts: &mut types::WaveFronts| {
        let wf_set: &mut Option<types::WaveFrontSet> = &mut wavefronts.wavefront_set[score];
        let wf_set: &mut types::WaveFrontSet = wf_set.as_mut().unwrap();

        let out_m_wf: &mut types::WaveFront = &mut wf_set.m.as_mut().unwrap();
        let out_d_wf: &mut types::WaveFront = &mut wf_set.d.as_mut().unwrap();

        let in_m_sub_wf: &types::WaveFront = awf_set.in_m_sub.unwrap();
        let in_m_gap_wf: &types::WaveFront = awf_set.in_m_gap.unwrap();
        let in_d_ext_wf: &types::WaveFront = awf_set.in_d_ext.unwrap();

        for k in lo..=hi {
            // Update D
            let out_k_index: usize = compute_k_index(out_m_wf.len(), out_m_wf.hi, k);
            let maybe_in_k_index_gap: Option<usize> = compute_maybe_k_index(in_m_gap_wf.len(), in_m_gap_wf.hi, k-1);
            let maybe_in_k_index_ext: Option<usize> = compute_maybe_k_index(in_d_ext_wf.len(), in_d_ext_wf.hi, k-1);
            let x: Option<i32> = match maybe_in_k_index_gap {
                Some(in_k_index_gap) => Some(in_m_gap_wf.offsets[in_k_index_gap]),
                None => None,
            };

            let y: Option<i32> = match maybe_in_k_index_ext {
                Some(in_k_index_ext) => Some(in_d_ext_wf.offsets[in_k_index_ext]),
                None => None,
            };

            let maybe_ins: Option<i32> = *vec![x, y].iter().max().unwrap();
            match maybe_ins {
                Some(v) => out_d_wf.offsets[out_k_index] = v+1,
                _ => out_d_wf.offsets[out_k_index] = -10,
            }


            // Update M
            let out_k_index = compute_k_index(out_m_wf.len(), out_m_wf.hi, k);
            let maybe_in_k_index = compute_maybe_k_index(in_m_sub_wf.len(), in_m_sub_wf.hi, k);

            match maybe_in_k_index {
                Some(in_k_index) => {
                    out_m_wf.offsets[out_k_index] = in_m_sub_wf.offsets[in_k_index] + 1;
                },
                None =>  {
                    out_m_wf.offsets[out_k_index] = -10;
                }
            };
        }
    };

    let assign_offsets_idm = |wavefronts: &mut types::WaveFronts| {
        let wf_set: &mut Option<types::WaveFrontSet> = &mut wavefronts.wavefront_set[score];
        let wf_set: &mut types::WaveFrontSet = wf_set.as_mut().unwrap();

        let out_m_wf: &mut types::WaveFront = &mut wf_set.m.as_mut().unwrap();
        let out_d_wf: &mut types::WaveFront = &mut wf_set.d.as_mut().unwrap();
        let out_i_wf: &mut types::WaveFront = &mut wf_set.i.as_mut().unwrap();

        eprintln!("awf_set: {:?}", awf_set);
        // let in_m_sub_wf: &types::WaveFront = awf_set.in_m_sub.unwrap();
        // let in_m_sub_wf: &types::WaveFront = awf_set.in_m_sub.unwrap();
        let maybe_in_m_sub_wf: Option<&types::WaveFront> = awf_set.in_m_sub;
        let maybe_in_m_gap_wf: Option<&types::WaveFront> = awf_set.in_m_gap;
        let maybe_in_d_ext_wf: Option<&types::WaveFront> = awf_set.in_d_ext;
        let maybe_in_i_ext_wf: Option<&types::WaveFront> = awf_set.in_i_ext;

        for k in lo..=hi {
            let out_k_index: usize = compute_k_index(out_m_wf.len(), out_m_wf.hi, k);

            // Update I

            if maybe_in_i_ext_wf.is_some() && maybe_in_m_gap_wf.is_some() {
                let in_i_ext_wf: &types::WaveFront = maybe_in_i_ext_wf.unwrap();
                let in_m_gap_wf: &types::WaveFront = maybe_in_m_gap_wf.unwrap();

                let maybe_in_k_index_gap: Option<usize> = compute_maybe_k_index(in_m_gap_wf.len(), in_m_gap_wf.hi, k-1);
                let maybe_in_k_index_ext: Option<usize> = compute_maybe_k_index(in_i_ext_wf.len(), in_i_ext_wf.hi, k-1);
                let x: Option<i32> = match maybe_in_k_index_gap {
                    Some(in_k_index_gap) => Some(in_m_gap_wf.offsets[in_k_index_gap]),
                    None => None,
                };

                let y: Option<i32> = match maybe_in_k_index_ext {
                    Some(in_k_index_ext) => Some(in_i_ext_wf.offsets[in_k_index_ext]),
                    None => None,
                };

                let maybe_ins: Option<i32> = *vec![x, y].iter().max().unwrap();
                match maybe_ins {
                    Some(v) => out_i_wf.offsets[out_k_index] = v+1,
                    _ => out_i_wf.offsets[out_k_index] = -10,
                }
            } else {
                out_i_wf.offsets[out_k_index] = -10;
            }


            // Update D
            if maybe_in_d_ext_wf.is_some() && maybe_in_m_gap_wf.is_some() {
                let in_d_ext_wf = maybe_in_d_ext_wf.unwrap();
                let in_m_gap_wf: &types::WaveFront = maybe_in_m_gap_wf.unwrap();

                let maybe_in_k_index_gap: Option<usize> = compute_maybe_k_index(in_m_gap_wf.len(), in_m_gap_wf.hi, k-1);
                let maybe_in_k_index_ext: Option<usize> = compute_maybe_k_index(in_d_ext_wf.len(), in_d_ext_wf.hi, k-1);
                let x: Option<i32> = match maybe_in_k_index_gap {
                    Some(in_k_index_gap) => Some(in_m_gap_wf.offsets[in_k_index_gap]),
                    None => None,
                };

                let y: Option<i32> = match maybe_in_k_index_ext {
                    Some(in_k_index_ext) => Some(in_d_ext_wf.offsets[in_k_index_ext]),
                    None => None,
                };

                let maybe_ins: Option<i32> = *vec![x, y].iter().max().unwrap();
                match maybe_ins {
                    Some(v) => out_d_wf.offsets[out_k_index] = v+1,
                    _ => out_d_wf.offsets[out_k_index] = -10,
                };

            } else {
                out_d_wf.offsets[out_k_index] = -10;
            }


            // Update M
            if maybe_in_m_sub_wf.is_some() {
                let in_m_sub_wf = maybe_in_m_sub_wf.unwrap();
                let out_k_index = compute_k_index(out_m_wf.len(), out_m_wf.hi, k);
                let maybe_in_k_index = compute_maybe_k_index(in_m_sub_wf.len(), in_m_sub_wf.hi, k);

                match maybe_in_k_index {
                    Some(in_k_index) => {
                        out_m_wf.offsets[out_k_index] = in_m_sub_wf.offsets[in_k_index] + 1;
                    },
                    None =>  {
                        out_m_wf.offsets[out_k_index] = -10;
                    }
                };
            } else {
                out_m_wf.offsets[out_k_index] = -10;
            }

        }
    };

    match wavefronts_to_allocate[..] {
        [types::WfType::M ] => { eprintln!("\t\tkernel: 0"); assign_offsets_m(wavefronts); },
        [types::WfType::M, types::WfType::I ] => { eprintln!("\t\tkernel: 2"); assign_offsets_im(wavefronts); },
        [types::WfType::M, types::WfType::D ] => { eprintln!("\t\tkernel: 1"); assign_offsets_dm(wavefronts); }
        [types::WfType::M, types::WfType::I, types::WfType::D ] => { eprintln!("\t\tkernel: 3"); assign_offsets_idm(wavefronts); },
        _ =>  { panic!("weird kernel: {:?}", wavefronts_to_allocate); }
    };
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
    assert_eq!(all_wavefronts.get_m_wavefront(score as i32).unwrap().get_offset(a_k as i32).cloned().unwrap(), 0);

    let mut match_posititons: Vec<(usize, usize, usize)> = Vec::new();

    // Print config
    if config.verbosity > 0 {
        eprintln!("Config\n\
                   \ttlen: {}\n\
                   \tqlen: {}\n\
                   \ta_k: {}\n\
                   \ta_offset: {}",
                  qlen, tlen, a_k, a_offset);
    }

    loop {
        // Extend the current wavefront
        if all_wavefronts.get_m_wavefront(score as i32).is_some() {
            let m_wf_mut: &mut types::WaveFront = &mut all_wavefronts
                .wavefront_set[score]
                .as_mut()
                .unwrap()
                .m
                .as_mut()
                .unwrap();
            wf_extend(m_wf_mut, match_lambda, &config, &mut match_posititons);
        }

        // Check whether we have reached the final point
        // Get the m-wavefront with the current score
        if utils::end_reached(all_wavefronts.get_m_wavefront(score as i32), a_k, a_offset) {

            /*
            print offsets
            ----

            eprintln!("k\tscore\toffset");
            for i in (0..=score).rev() {
                let m_wf = all_wavefronts.get_m_wavefront(i as i32);

                match m_wf {
                    Some (wf) => {
                        for k in wf.lo..=wf.hi {
                            let offset = wf.get_offset(k);
                            eprintln!("{}\t{}\t{:?}", k, i, offset);
                        }
                    }
                    _ => {}
                };
                eprintln!();

            }
             */
            /*
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
             */
            break
        }

        score += 1;

        // compute the next wavefront
        wf_next(&mut all_wavefronts, score, config);
    }

    //let cigar = String::new();
    let cigar = wf_traceback(&all_wavefronts, score, config);

    let each_wf = vec![ types::WfType::M ];
    utils::debug_utils::visualize_all(&all_wavefronts, a_offset, &each_wf, &match_posititons, config, score);

    (score, cigar)
}


#[cfg(test)]
mod tests {
    mod test_config {
        pub static CONFIG: crate::types::Config = crate::types::Config {
            adapt: false,
            verbosity: 5,
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
                let text  = "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";
                let query = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT";

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

                let (score, cigar) = wf_align(tlen as u32, qlen as u32, &match_lambda, &test_config::CONFIG);
                // eprintln!("Result:\n\tScore: {} Cigar {}", score, cigar);
                // crate::utils::backtrace_utils::print_aln(&cigar[..], t, q);
            }
        }
    }
}
