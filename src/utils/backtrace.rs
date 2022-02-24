// TODO: remove
use super::*;
// O(n)
pub fn print_aln(cigar: &str, t: &[u8], q: &[u8]) {
    let mut query = String::new();
    let mut marker = String::new();
    let mut text = String::new();

    let space = &String::from(" ")[..];
    let dash = &String::from("-")[..];
    let vertical_bar = &String::from("|")[..];

    let mut q_iter = q.iter();
    let mut t_iter = t.iter();

    // TODO: does this unwrap or gap introduce a bug?
    let foo = |y: Option<&u8>| -> String { format!("{}", *y.unwrap_or(&b'-') as char) };

    for c in cigar.as_bytes() {
        match c {
            b'M' => {
                query.push_str(&foo(q_iter.next()));
                marker.push_str(vertical_bar);
                text.push_str(&foo(t_iter.next()));
            }

            b'X' => {
                query.push_str(&foo(q_iter.next()));
                marker.push_str(space);
                text.push_str(&foo(t_iter.next()));
            }

            b'I' => {
                query.push_str(dash);
                marker.push_str(space);
                text.push_str(&foo(t_iter.next()));
            }

            b'D' => {
                query.push_str(&foo(q_iter.next()));
                marker.push_str(space);
                text.push_str(dash);
            }

            _ => {
                panic!("[utils::backtrace_utils::print_aln] found char not M, I, X or D")
            }
        }
    }

    eprintln!();
    eprintln!("{}", query);
    eprintln!("{}", marker);
    eprintln!("{}", text);
}

// Compare current to next and accumulate counts
// O(n)
pub fn run_length_encode(cigar: &str, reverse: bool) -> String {
    let mut cigar = String::from(cigar);
    if reverse {
        cigar = cigar.chars().rev().collect::<String>();
    }

    let mut xcigar = String::new();

    // edge cases
    if cigar.len() == 0 {
        panic!("[wfa::utils::run_length_encode] empty cigar");
    } else if cigar.len() == 1 {
        xcigar.push_str(&format!("{}{}", 1, cigar));
        xcigar
    } else {
        let mut chars = cigar.chars();
        let mut current: char = chars.next().unwrap();
        let mut count = 1;

        for c in chars {
            if c == current {
                count += 1;
            } else {
                xcigar.push_str(&format!("{}{}", count, current));
                current = c;
                count = 1;
            }
        }

        // last run
        xcigar.push_str(&format!("{}{}", count, current));

        xcigar
    }
}

pub fn wflambda_backtrace_matches_check<G>(
    offset: &mut i32,
    cigar: &mut String,
    num_matches: u32,
    k: i32,
    traceback_lambda: &mut G
) where
    G: FnMut((i32, i32), (i32, i32)) {

    {
        // let o = *offset as u64;
        let query_stop = compute_v(*offset, k);
        let target_stop = compute_h(*offset, k);


        let query_start = query_stop - num_matches as i32;
        let target_start = target_stop - num_matches as i32;

        let query = (query_start as i32, query_stop as i32);
        let target = (target_start as i32, target_stop as i32);

        traceback_lambda(query, target);
    }

    // TODO: improve this add M x-times and subtruct offset by num_matches
    (0..num_matches).for_each(|_| {
        // let v = compute_v(*offset, k, central_diagonal);
        // let h = compute_h(*offset, k, central_diagonal);

        cigar.push('M');
        *offset -= 1;
    });
}

// TODO: will this ever run in regions without a match?
pub fn backtrace_matches_check(offset: &mut i32, cigar: &mut String, num_matches: u32) {
    cigar.extend(repeat_char('M', num_matches));
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
