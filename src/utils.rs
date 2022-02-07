use num;
use ndarray::{Array2, Array};
use ndarray_to_img;
use std::fs;
use std::path;

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


pub fn compute_v(offset: i32, k: i32) -> i32 {
    offset - k
}

pub fn compute_h(offset: i32, _: i32) -> i32 {
    offset
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

        let foo = |y: Option<&u8>| -> String {
            format!("{}", *y.unwrap_or( &b'X' ) as char)
        };

        for c in cigar.as_bytes() {
            match c {
                b'M'=> {
                    query.push_str( &foo(q_iter.next()) );
                    marker.push_str(vertical_bar);
                    text.push_str( &foo(t_iter.next()) );
                },

                b'X' => {
                    query.push_str( &foo(q_iter.next()) );
                    marker.push_str(space);
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

    pub fn backtrace_matches_check(
        offset: &mut i32,
        cigar: &mut String,
        num_matches: u32,
        k: i32,
    ) {
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
        score: i32,
        k: i32
    ) -> Option<i32> {
        if score < 0 {
            return None;
        }

        // if m_wf.lo <= k + 1 && k + 1 <= m_wf.hi {
        all_wavefronts
            .get_m_wavefront(score)
            .and_then(|d_wf| d_wf.get_offset(k+1))
            .cloned()
    }

    pub fn backtrace_deletion_extend_offset(
        all_wavefronts: &crate::types::WaveFronts,
        score: i32,
        k: i32
    ) -> Option<i32> {
        if score < 0 {
            return None;
        }

        // is d_wf.lo <= k + 1 && k + 1 <= d_wf.hi
        all_wavefronts
            .get_d_wavefront(score)
            .and_then(|d_wf| d_wf.get_offset(k+1))
            .cloned()
    }

    pub fn backtrace_insertion_open_offset(
        all_wavefronts: &crate::types::WaveFronts,
        score: i32,
        k: i32
    ) -> Option<i32> {
        if score < 0 {
            return None;
        }

        all_wavefronts
            .get_m_wavefront(score)
            .and_then(|i_wf| i_wf.get_offset(k-1))
            .cloned()
            .map(|x| x + 1)
    }

    pub fn backtrace_insertion_extend_offset(
        all_wavefronts: &crate::types::WaveFronts,
        score: i32,
        k: i32
    ) -> Option<i32> {
        if score < 0 {
            return None;
        }

        // if i_wf.lo <= k - 1 && k - 1 <= i_wf.hi {
        all_wavefronts
            .get_i_wavefront(score)
            .and_then(|i_wf| i_wf.get_offset(k-1))
            .cloned()
            .map(|x| x + 1)
    }

    pub fn backtrace_mismatch_offset(
        all_wavefronts: &crate::types::WaveFronts,
        score: i32,
        k: i32
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
}

pub mod debug_utils {
    use super::*;


    pub fn visualize_all(
        all_wavefronts: &types::WaveFronts,
        a_offset: u32,
        each_wf: &Vec<types::WfType>,
        match_positions: &Vec<(i32, i32, usize)>,
        config: &types::Config,
        score: usize
    ) {
        if config.verbosity > 1 {
            eprintln!("[utils::visualize_all]");
        }

        if config.verbosity > 2 {
            eprintln!("\tPopulating matrix");
        }


        let dim = (a_offset as usize+1, a_offset as usize+1);
        let x = ndarray::Dim(dim);
        let mut matrix: Array2<Option<i32>> = Array::from_elem(x, None);

        eprintln!("\t\tk\tscore\toffset\t(v,h)");
        for s in (0..=score).rev() {
            let wf_specific: &types::WaveFront = match all_wavefronts.get_m_wavefront(s as i32) {
                Some(m) => m,
                _ => continue,
            };

            let lo = wf_specific.lo;
            let hi = wf_specific.hi;
            let offsets = &wf_specific.offsets;
            let len = offsets.len();

            // eprintln!("\t\tscore: {} lo {} hi {}", s, lo, hi);

            for k in lo..=hi {
                let k_index: usize = compute_k_index(len, k, hi);
                let m_s_k: i32 = offsets[k_index];

                for offset in  m_s_k..=m_s_k {
                    let v = compute_v(offset, k);
                    let h = compute_h(offset, k);

                    // eprintln!("offset: {}\tk: {}\tscore: {}\t({}, {})", m_s_k,  k, score, v, h);

                    eprintln!("\t\t{}\t{}\t{}\t({},{})", k, s, offset, v, h);


                    if v < 0 || h < 0 || v >= dim.0 as i32 || h >= dim.0 as i32 {
                        continue;
                    }

                    let v = v as usize;
                    let h = h as usize;

                    if config.verbosity > 5 {
                        // eprintln!("\t\t({},{})\t{}\t{}", v, h, s, offset);
                    }


                    matrix[[v,h]] = Some(offset);
                }
            }

            eprintln!();
        }

        for (v, h, score) in match_positions {
            let v = *v;
            let h = *h;

            if v < 0 || h < 0 || v >= dim.0 as i32 || h >= dim.0 as i32 {
                continue;
            }

            let v = v as usize;
            let h = h as usize;

            if v == h {
                matrix[[v, h]] = Some(*score as i32);
            }
        }

        gen_image(&matrix, config);
    }

    fn gen_image(matrix: &Array2<Option<i32>>, config: &types::Config) {
        if config.verbosity > 1 {
            eprintln!("[utils::gen_image]");
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

        let out_dir = path::Path::new("./debug/out/");
        // assume this will never fail
        if !out_dir.exists() {
            fs::create_dir_all(out_dir).unwrap();
        }

        let image_path_buf: path::PathBuf = out_dir.join("all.png");
        let image_name: &str = image_path_buf.to_str().unwrap();
        ndarray_to_img::generate_image(&scaled_matrix, &config, image_name).unwrap();
    }
}
