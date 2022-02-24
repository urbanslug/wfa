use crate::types;

fn compute_wf_next_limits(
    wavefronts: &types::WaveFronts,
    score: usize,
    config: &types::Config,
) -> (Option<i32>, Option<i32>) {
    enum WfLimit {
        Hi,
        Lo,
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
                WfLimit::Lo => Some(wf.lo),
            },
            _ => None,
        }
    };

    // what if a hi is empty?
    let hi: Option<i32> = vec![
        foo(wavefronts.get_m_wavefront(s_x), WfLimit::Hi),
        foo(wavefronts.get_m_wavefront(s_o_e), WfLimit::Hi),
        foo(wavefronts.get_i_wavefront(s_e), WfLimit::Hi),
        foo(wavefronts.get_d_wavefront(s_e), WfLimit::Hi),
    ]
    .iter()
    .max()
    .unwrap()
    .map(|x| x + 1);

    let maybe_los: Vec<Option<i32>> = vec![
        foo(wavefronts.get_m_wavefront(s_x), WfLimit::Lo),
        foo(wavefronts.get_m_wavefront(s_o_e), WfLimit::Lo),
        foo(wavefronts.get_i_wavefront(s_e), WfLimit::Lo),
        foo(wavefronts.get_d_wavefront(s_e), WfLimit::Lo),
    ]
    .iter()
    .filter(|x| x.is_some())
    .cloned()
    .collect();

    if maybe_los.is_empty() {
        return (hi, None);
    }

    let lo: Option<i32> = maybe_los.iter().min().unwrap().map(|x| x - 1);

    (hi, lo)
}

// TODO: move to types
#[derive(Debug)]
struct AWFSet<'a> {
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
    config: &'a types::Config,
) -> AWFSet<'a> {
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
    if awf_set.in_m_gap.is_some() || awf_set.in_i_ext.is_some() {
        wavefronts_to_allocate.push(types::WfType::I);
    }

    // Allocate D-Wavefront
    if awf_set.in_m_gap.is_some() || awf_set.in_d_ext.is_some() {
        wavefronts_to_allocate.push(types::WfType::D);
    }

    wavefronts_to_allocate
}

fn foobar<'a>(
    wavefronts: &'a mut types::WaveFronts,
    awf_set: &AWFSet,
    lo: i32,
    hi: i32,
    score: usize,
) -> Vec<types::WfType> {
    let mut wavefronts_to_allocate = vec![types::WfType::M];

    let maybe_out_m_wf = Some(types::WaveFront::new(hi, lo));
    let mut maybe_out_i_wf = None;
    let mut maybe_out_d_wf = None;

    // Allocate I-Wavefront
    if awf_set.in_m_gap.is_some() || awf_set.in_i_ext.is_some() {
        maybe_out_i_wf = Some(types::WaveFront::new(hi, lo));
        wavefronts_to_allocate.push(types::WfType::I);
    }

    // Allocate D-Wavefront
    if awf_set.in_m_gap.is_some() || awf_set.in_d_ext.is_some() {
        maybe_out_d_wf = Some(types::WaveFront::new(hi, lo));
        wavefronts_to_allocate.push(types::WfType::D);
    }

    let prev_score = wavefronts.wavefront_set.len();

    for s in prev_score..=score {
        if s == score {
            let wf = types::WaveFrontSet {
                i: maybe_out_i_wf.clone(),
                m: maybe_out_m_wf.clone(),
                d: maybe_out_d_wf.clone(),
            };

            wavefronts.wavefront_set.push(Some(wf));
        } else {
            wavefronts.wavefront_set.push(None);
        }
    }

    // awf_set.out_m = wavefronts.get_m_wavefront(score as i32);

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
    let mut awf_set = fetch_wf(score, &cloned_wf, config);

    if awf_set.in_m_sub.is_none()
        && awf_set.in_m_gap.is_none()
        && awf_set.in_d_ext.is_none()
        && awf_set.in_d_ext.is_none()
    {
        if config.verbosity > 4 {
            eprintln!("\t\tskipping score {}", score);
        }
        return;
    }

    if config.verbosity > 5 {
        eprintln!(
            "\t\ts {} s - o - e {} s - e {} s - x {}",
            s, unsigned_s_o_e, unsigned_s_e, unsigned_s_x
        );
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

    // eprintln!("\t\tscore: {} lo {} hi {}", score, lo, hi);

    // Allocate the next wave front
    let wavefronts_to_allocate = alloc_wf(&awf_set);
    // alloc_wf(&wavefronts, &awf_set, score, lo, hi);
    // let allocation_result = wavefronts.allocate_wavefronts(score as u32, lo, hi, &wavefronts_to_allocate);
    /*
    match allocation_result {
        Ok(_) => {
            if config.verbosity > 254 {
                eprintln!("\t\tallocated wfs for score: {} \n\t {:?}", score, wavefronts.wavefront_set[score])
            }
        },
        Err(msg) => { panic!("{}", msg) }
    };
    */

    let wavefronts_to_allocate = foobar(wavefronts, &awf_set, lo, hi, score);
    // awf_set.out_m = wavefronts.get_m_wavefront(score as i32);
    // awf_set.out_i = wavefronts.get_i_wavefront(score as i32);
    // awf_set.out_d = wavefronts.get_d_wavefront(score as i32);

    // let wave_length: usize = utils::compute_wave_length(lo, hi);

    let compute_k_index = |wavefront_length: usize, hi: i32, k: i32| -> usize {
        wavefront_length - ((hi - k) as usize) - 1
    };

    let compute_maybe_k_index = |wavefront_length: usize, hi: i32, k: i32| -> Option<usize> {
        if hi < k {
            return None;
        };

        let candidate_index = wavefront_length as i32 - (hi - k) - 1;

        if candidate_index < 0 {
            return None;
        }

        Some(candidate_index as usize)
    };

    let affine_wavefront_cond_fetch = |wf: &types::WaveFront, k: i32| -> i32 {
        // eprintln!("{:?}", wf);
        compute_maybe_k_index(wf.len(), wf.hi, k)
            .map(|k_index: usize| wf.offsets[k_index])
            .unwrap_or(-10)
    };

    let assign_offsets_m = |wavefronts: &mut types::WaveFronts| {
        let wf_set: &mut Option<types::WaveFrontSet> = &mut wavefronts.wavefront_set[score];
        let wf_set: &mut types::WaveFrontSet = wf_set.as_mut().unwrap();
        let out_m_wf: &mut types::WaveFront = &mut wf_set.m.as_mut().unwrap();
        let wavefront_len = out_m_wf.len();

        let in_m_wf: &types::WaveFront = awf_set.in_m_sub.unwrap();

        for k in lo..=hi {
            let k_index = compute_k_index(wavefront_len, out_m_wf.hi, k);
            let mut offset = affine_wavefront_cond_fetch(in_m_wf, k);
            if offset != -10 {
                offset += 1
            }
            out_m_wf.offsets[k_index] = offset;
            // eprintln!("\t\tk {}\tM {}", k, offset);
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
            // comapre gap open on M and gap extend on I
            let k_index: usize = compute_k_index(out_i_wf.len(), out_i_wf.hi, k);

            let x: i32 = affine_wavefront_cond_fetch(in_m_gap_wf, k - 1);
            let y: i32 = affine_wavefront_cond_fetch(in_i_ext_wf, k - 1);
            let ins: i32 = *vec![x, y].iter().max().unwrap();

            out_i_wf.offsets[k_index] = ins + 1;

            // Update M
            let k_index = compute_k_index(out_m_wf.len(), out_m_wf.hi, k);
            let sub: i32 = affine_wavefront_cond_fetch(in_m_sub_wf, k) + 1;
            let sub: i32 = *vec![sub, ins].iter().max().unwrap();

            out_m_wf.offsets[k_index] = sub;
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
            // comapre gap open on M and gap extend on I
            let k_index: usize = compute_k_index(out_d_wf.len(), out_d_wf.hi, k);

            let x: i32 = affine_wavefront_cond_fetch(in_m_gap_wf, k - 1);
            let y: i32 = affine_wavefront_cond_fetch(in_d_ext_wf, k - 1);
            let del = *vec![x, y].iter().max().unwrap();

            out_d_wf.offsets[k_index] = del;

            // Update M
            let k_index = compute_k_index(out_m_wf.len(), out_m_wf.hi, k);
            let sub: i32 = affine_wavefront_cond_fetch(in_m_sub_wf, k) + 1;
            let max_m: i32 = *vec![sub, del].iter().max().unwrap();

            out_m_wf.offsets[k_index] = sub;
        }
    };

    let assign_offsets_idm = |wavefronts: &mut types::WaveFronts| {
        let wf_set: &mut Option<types::WaveFrontSet> = &mut wavefronts.wavefront_set[score];
        let wf_set: &mut types::WaveFrontSet = wf_set.as_mut().unwrap();

        let out_m_wf: &mut types::WaveFront = &mut wf_set.m.as_mut().unwrap();
        let out_d_wf: &mut types::WaveFront = &mut wf_set.d.as_mut().unwrap();
        let out_i_wf: &mut types::WaveFront = &mut wf_set.i.as_mut().unwrap();

        // eprintln!("{:#?}", awf_set);

        let maybe_in_m_sub_wf: Option<&types::WaveFront> = awf_set.in_m_sub;
        let maybe_in_m_gap_wf: Option<&types::WaveFront> = awf_set.in_m_gap;
        let maybe_in_d_ext_wf: Option<&types::WaveFront> = awf_set.in_d_ext;
        let maybe_in_i_ext_wf: Option<&types::WaveFront> = awf_set.in_i_ext;

        // compute min_hi
        let min_hi: Option<i32> = vec![
            maybe_in_d_ext_wf,
            maybe_in_i_ext_wf,
            maybe_in_m_gap_wf,
            maybe_in_m_sub_wf,
        ]
        .into_iter()
        .filter(|x| x.is_some())
        .map(|maybe_wf| maybe_wf.map(|wf| wf.hi))
        .min()
        .unwrap();
        let min_hi = min_hi.unwrap();

        let max_lo: Option<i32> = vec![
            maybe_in_d_ext_wf,
            maybe_in_i_ext_wf,
            maybe_in_m_gap_wf,
            maybe_in_m_sub_wf,
        ]
        .into_iter()
        .filter(|x| x.is_some())
        .map(|maybe_wf| maybe_wf.map(|wf| wf.lo))
        .max()
        .unwrap();

        let max_lo = max_lo.unwrap();

        for k in lo..=hi {
            // Update I
            let k_index: usize = compute_k_index(out_i_wf.len(), out_i_wf.hi, k);
            let ins_m =
                maybe_in_m_gap_wf.and_then(|m_gap| Some(affine_wavefront_cond_fetch(m_gap, k - 1)));
            let ins_i =
                maybe_in_i_ext_wf.and_then(|i_ext| Some(affine_wavefront_cond_fetch(i_ext, k - 1)));
            let ins: i32 = vec![ins_m, ins_i]
                .into_iter()
                .max()
                .unwrap()
                .map(|i| i + 1)
                .unwrap();
            // let ins: i32 = maybe_ins.unwrap_or(-10);
            out_i_wf.offsets[k_index] = ins;

            // Update D
            let k_index: usize = compute_k_index(out_d_wf.len(), out_d_wf.hi, k);
            let del_m =
                maybe_in_m_gap_wf.and_then(|m_gap| Some(affine_wavefront_cond_fetch(m_gap, k + 1)));
            let del_i =
                maybe_in_d_ext_wf.and_then(|d_ext| Some(affine_wavefront_cond_fetch(d_ext, k + 1)));
            let del: i32 = vec![del_m, del_i].into_iter().max().unwrap().unwrap();
            // let del: i32 = maybe_del.unwrap_or(-10);
            out_d_wf.offsets[k_index] = del;

            // Update M
            let k_index: usize = compute_k_index(out_m_wf.len(), out_m_wf.hi, k);
            let sub_m: Option<i32> = maybe_in_m_sub_wf
                .and_then(|m_sub| Some(affine_wavefront_cond_fetch(m_sub, k)))
                .map(|x| if x == -10 { -10 } else { x + 1 });
            let sub: i32 = vec![sub_m, Some(ins), Some(del)]
                .into_iter()
                .max()
                .unwrap()
                .unwrap();
            // let sub = maybe_sub.unwrap_or(-10);
            out_m_wf.offsets[k_index] = sub;

            // eprintln!("\t\tk {}\tM {}\tI {}\tD {}", k, sub, ins, del);
        }
        /*
        for k in max_lo..=min_hi {

            // Update I
            // comapre gap open on M and gap extend on I
            let k_index: usize = compute_k_index(out_i_wf.len(), out_i_wf.hi, k);

            let x: Option<i32> = affine_wavefront_cond_fetch(in_m_gap_wf, k-1);
            let y: Option<i32> = affine_wavefront_cond_fetch(in_i_ext_wf, k-1);
            let maybe_ins: Option<i32> = *vec![x, y].iter().max().unwrap();
            let ins = match maybe_ins {
                Some(v) => v+1,
                None => -10,
            };

            out_i_wf.offsets[k_index] = ins;

            // Update D
            // comapre gap open on M and gap extend on I
            let k_index: usize = compute_k_index(out_d_wf.len(), out_d_wf.hi, k);

            let x: Option<i32> = affine_wavefront_cond_fetch(in_m_gap_wf, k-1);
            let y: Option<i32> = affine_wavefront_cond_fetch(in_d_ext_wf, k-1);
            let maybe_del: Option<i32> = *vec![x, y].iter().max().unwrap();
            let del = match maybe_del {
                Some(v) => v,
                None => -10,
            };

            out_d_wf.offsets[k_index] = del;

            // Update M
            let k_index = compute_k_index(out_m_wf.len(), out_m_wf.hi, k);
            let maybe_sub: Option<i32> = affine_wavefront_cond_fetch(in_m_sub_wf, k);
            let maybe_sub: Option<i32> = *vec![maybe_sub, Some(del), Some(ins)].iter().max().unwrap();

            let sub = match maybe_sub {
                Some(v) => v,
                None => -10,
            };

            out_m_wf.offsets[k_index] = sub;
        }
        for k in min_hi..hi {}
        */
    };

    match wavefronts_to_allocate[..] {
        [types::WfType::M] => {
            // eprintln!("\t\tkernel: 0");
            assign_offsets_m(wavefronts);
        }
        [types::WfType::M, types::WfType::I] => {
            // eprintln!("\t\tkernel: 2");
            assign_offsets_im(wavefronts);
        }
        [types::WfType::M, types::WfType::D] => {
            // eprintln!("\t\tkernel: 1");
            assign_offsets_dm(wavefronts);
        }
        [types::WfType::M, types::WfType::I, types::WfType::D] => {
            // eprintln!("\t\tkernel: 3");
            assign_offsets_idm(wavefronts);
        }
        _ => {
            panic!("weird kernel: {:?}", wavefronts_to_allocate);
        }
    };
}

#[cfg(test)]
mod tests {
}
