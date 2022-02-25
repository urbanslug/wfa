use crate::utils;

pub fn wflambda_backtrace_matches_check<G>(
    offset: &mut i32,
    cigar: &mut String,
    num_matches: u32,
    k: i32,
    traceback_lambda: &mut G
) where
    G: FnMut((i32, i32), (i32, i32)) -> bool,
{
    let query_stop = utils::compute_v(*offset, k);
    let target_stop = utils::compute_h(*offset, k);

    let query_start = query_stop - num_matches as i32;
    let target_start = target_stop - num_matches as i32;

    let query = (query_start as i32, query_stop as i32);
    let target = (target_start as i32, target_stop as i32);

    if traceback_lambda(query, target) {
        cigar.extend(utils::repeat_char('M', num_matches));
        *offset -= num_matches as i32;
    }
}

// TODO: will this ever run in regions without a match?
pub fn backtrace_matches_check(
    offset: &mut i32,
    cigar: &mut String,
    num_matches: u32
) {
    cigar.extend(utils::repeat_char('M', num_matches));
    *offset -= num_matches as i32;
}

pub fn backtrace_deletion_open_offset(
    all_wavefronts: &crate::types::WaveFronts,
    score: i32,
    k: i32,
) -> Option<i32> {
    if score < 0 {
        return None;
    }

    // if m_wf.lo <= k + 1 && k + 1 <= m_wf.hi {
    all_wavefronts
        .get_m_wavefront(score)
        .and_then(|d_wf| d_wf.get_offset(k + 1))
        .cloned()
}

pub fn backtrace_deletion_extend_offset(
    all_wavefronts: &crate::types::WaveFronts,
    score: i32,
    k: i32,
) -> Option<i32> {
    if score < 0 {
        return None;
    }

    // is d_wf.lo <= k + 1 && k + 1 <= d_wf.hi
    all_wavefronts
        .get_d_wavefront(score)
        .and_then(|d_wf| d_wf.get_offset(k + 1))
        .cloned()
}

pub fn backtrace_insertion_open_offset(
    all_wavefronts: &crate::types::WaveFronts,
    score: i32,
    k: i32,
) -> Option<i32> {
    if score < 0 {
        return None;
    }

    all_wavefronts
        .get_m_wavefront(score)
        .and_then(|i_wf| i_wf.get_offset(k - 1))
        .cloned()
        .map(|x| x + 1)
}

pub fn backtrace_insertion_extend_offset(
    all_wavefronts: &crate::types::WaveFronts,
    score: i32,
    k: i32,
) -> Option<i32> {
    if score < 0 {
        return None;
    }

    // if i_wf.lo <= k - 1 && k - 1 <= i_wf.hi {
    all_wavefronts
        .get_i_wavefront(score)
        .and_then(|i_wf| i_wf.get_offset(k - 1))
        .cloned()
        .map(|x| x + 1)
}

pub fn backtrace_mismatch_offset(
    all_wavefronts: &crate::types::WaveFronts,
    score: i32,
    k: i32,
) -> Option<i32> {
    if score < 0 {
        return None;
    }

    all_wavefronts
        .get_m_wavefront(score)
        .and_then(|m_wf| m_wf.get_offset(k))
        .cloned()
        .map(|x| x + 1)
}
