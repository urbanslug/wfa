use num;
use ndarray::{Array2, Array};
use ndarray_to_img;

use crate::types;



// How many cells does the wave cross?
pub fn compute_wave_length(lo: i32, hi: i32) -> usize {
    (hi-lo+1) as usize
}

pub fn to_usize_or_zero<T: num::cast::ToPrimitive>(n: T) -> usize {
    n.to_usize().unwrap_or(0)
}

pub fn compute_k_index(length: usize, k: i32, hi: i32) -> usize {
    // we expect hi - k to always be +ve
    length - ((hi - k) as usize) - 1
}

pub fn sub_else_zero(lhs: isize, rhs: isize) -> isize {
    let result: isize = lhs - rhs;
    if result < 0 {
        0
    } else {
        result
    }
}

pub fn abs_sub(lhs: i32, rhs: i32) -> i32 {
    let result = lhs - rhs;
    num::abs(result)
}

pub fn compute_v(offset: i32, k: i32) -> usize {
    abs_sub(offset, k as i32) as usize
}

pub fn compute_h(offset: i32, _: i32) -> usize {
    offset as usize
}

pub fn end_reached(
    m_wavefront: Option<&types::WaveFront>,
    a_k: usize,
    a_offset: u32
) -> bool {
    let m_wavefront = match m_wavefront {
        Some(wf) => wf,
        _ => { return false }
    };

    let k_index: usize = compute_k_index(m_wavefront.len(), a_k as i32, m_wavefront.hi);
    let m_s_k = m_wavefront.offsets[k_index]; // m_offset

    m_s_k > 0 && (m_s_k as u32) >= a_offset
}

pub mod backtrace_utils {

    pub fn print_aln(cigar: &str, t: &[u8], q: &[u8]) {

        let mut query = String::new();
        let mut marker = String::new();
        let mut text = String::new();

        let space = &String::from(" ")[..];
        let dash = &String::from("-")[..];
        let vertical_bar = &String::from("|")[..];

        let mut q_iter = q.iter();
        let mut t_iter = t.iter();

        let foo = |y: Option<&u8>| -> String {
            format!("{}", *y.unwrap_or( &b'X' ) as char)
        };

        for c in cigar.as_bytes() {
            match c {
                b'M' | b'X' => {
                    query.push_str( &foo(q_iter.next()) );
                    marker.push_str(vertical_bar);
                    text.push_str( &foo(t_iter.next()) );
                },

                b'I' => {
                    query.push_str(dash);
                    marker.push_str(space);
                    text.push_str( &foo(t_iter.next()) );
                },

                b'D' => {
                    query.push_str( &foo(q_iter.next()) );
                    marker.push_str(space);
                    text.push_str(dash);
                },

                _ => { panic!("[utils::backtrace_utils::print_aln] found char not M, I, X or D") },
            }
        }

        eprintln!("{}", cigar);

        eprintln!();

        eprintln!("{}", query);
        eprintln!("{}", marker);
        eprintln!("{}", text);
    }

    // just make the cigar proper
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

    pub fn backtrace_matches_check(
        offset: &mut isize,
        cigar: &mut String,
        num_matches: usize,
        k: i32,
    )
    {
        // TODO: improve this add M x-times and subtruct offset by num_matches
        (0..num_matches).for_each(|_| {
            // let v = compute_v(*offset, k, central_diagonal);
            // let h = compute_h(*offset, k, central_diagonal);

            cigar.push('M');
            *offset -= 1;
        });
    }

    pub fn backtrace_deletion_open_offset(
        all_wavefronts: &crate::types::WaveFronts,
        score: isize,
        k: i32
    ) -> Option<isize> {
        if score < 0 {
            return None;
        }

        // if m_wf.lo <= k + 1 && k + 1 <= m_wf.hi {
        all_wavefronts
            .get_d_wavefront(score as usize)
            .and_then(|d_wf| d_wf.get_offset(k+1))
            .cloned()
            .map(|x| x as isize)
    }

    pub fn backtrace_deletion_extend_offset(
        all_wavefronts: &crate::types::WaveFronts,
        score: isize,
        k: i32
    ) -> Option<isize> {
        if score < 0 {
            return None;
        }

        // is d_wf.lo <= k + 1 && k + 1 <= d_wf.hi
        all_wavefronts
            .get_d_wavefront(score as usize)
            .and_then(|d_wf| d_wf.get_offset(k+1))
            .cloned()
            .map(|x| x as isize)

    }

    pub fn backtrace_insertion_open_offset(
        all_wavefronts: &crate::types::WaveFronts,
        score: isize,
        k: i32
    ) -> Option<isize> {
        if score < 0 {
            return None;
        }

        all_wavefronts
            .get_i_wavefront(score as usize)
            .and_then(|i_wf| i_wf.get_offset(k-1))
            .cloned()
            .map(|x| x as isize + 1)
    }

    pub fn backtrace_insertion_extend_offset(
        all_wavefronts: &crate::types::WaveFronts,
        score: isize,
        k: i32
    ) -> Option<isize> {
        if score < 0 {
            return None;
        }

        // if i_wf.lo <= k - 1 && k - 1 <= i_wf.hi {
        all_wavefronts
            .get_i_wavefront(score as usize)
            .and_then(|i_wf| i_wf.get_offset(k-1))
            .cloned()
            .map(|x| x as isize + 1)
    }

    pub fn backtrace_mismatch_offset(
        all_wavefronts: &crate::types::WaveFronts,
        score: isize,
        k: i32
    ) -> Option<isize> {
        if score < 0 {
            return None;
        }

        all_wavefronts
            .get_m_wavefront(score as usize)
            .and_then(|m_wf| m_wf.get_offset(k))
            .cloned()
            .map(|x| x as isize + 1)
    }
}

pub mod debug_utils {
    use super::*;

    pub fn visualize_all(
        all_wavefronts: &types::WaveFronts,
        a_offset: u32,
        each_wf: &Vec<types::WfType>
    ) {
        let dim = (a_offset as usize+10, a_offset as usize+10);
        let x = ndarray::Dim(dim);
        let mut matrix: Array2<Option<i32>> = Array::from_elem(x, None);

        // Loop through all the wavefront sets
        for wf_set in all_wavefronts.wavefront_set.iter() {

            let wf_set = match wf_set {
                Some(w) => w,
                None => { continue; },
            };

            // if there's something in here
            for wf in each_wf {
                let wf_specific = match wf {
                    types::WfType::I => { match &wf_set.i {
                        Some(w) => w,
                        _ => { continue; }
                    }},
                    types::WfType::D => { match &wf_set.d {
                        Some(w) => w,
                        _ => { continue; }
                    }},
                    types::WfType::M => { match &wf_set.m {
                        Some(w) => w,
                        _ => { continue; }
                    }},
                };

                let lo = wf_specific.lo;
                let hi = wf_specific.hi;
                let offsets = &wf_specific.offsets;
                let len = offsets.len();

                for k in lo..=hi {

                    let k_index: usize = compute_k_index(len, k, hi);
                    let m_s_k: i32 = offsets[k_index];

                    for offset in  0..=m_s_k {
                        let v: usize = compute_v(offset, k);
                        let h: usize = compute_h(offset, k);

                        // eprintln!("offset: {}\tk: {}\tscore: {}\t({}, {})", m_s_k,  k, score, v, h);

                        if v >= dim.0 || h >= dim.0 {
                            // continue;
                        }


                        matrix[[v,h]] = Some(offset);
                    }
                }
            }
        }

        let config = ndarray_to_img::Config {
            verbosity: 0,
            with_color: true,
            annotate_image: true,
            draw_diagonal: true,
            draw_boundaries: true,
            scaling_factor: 10,
        };

        let scaled_matrix = ndarray_to_img::scale_matrix(&matrix, &config);
        let image_name = "all.png";
        ndarray_to_img::generate_image(&scaled_matrix, &config, &image_name).unwrap();
    }
}

