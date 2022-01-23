use num;
use ndarray::{Array2};
use crate::types;

use ndarray_to_img;

// How many cells does the wave cross?
pub fn compute_wave_length(lo: i32, hi: i32) -> usize {
		(hi-lo+1) as usize
}

// 
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

		pub fn backtrace_deletion_extend_offset(
				all_wavefronts: &crate::types::WaveFronts,
				score: isize,
				k: i32
		) -> Option<isize> {
				if score < 0 {
						return None;
				}

				// we know it's >= 0
				let score = score as usize;

				// make safe use get to return None
				let d_wf = &all_wavefronts.wavefront_set[score].d;

				if d_wf.lo <= k + 1 && k + 1 <= d_wf.hi {
            let k_index = crate::utils::compute_k_index(d_wf.vals.len(), k, d_wf.hi);
						// make safe
						let d_s_k = d_wf.vals[k_index+1];
						Some(d_s_k as isize)
        } else {
            None
        }
		}

		pub fn backtrace_deletion_open_offset(
				all_wavefronts: &crate::types::WaveFronts,
				score: isize,
				k: i32
		) -> Option<isize> {
				if score < 0 {
						return None;
				}

				// we know it's >= 0
				let score = score as usize;

				// make safe use get to return None
				let m_wf = &all_wavefronts.wavefront_set[score].m;

				if m_wf.lo <= k + 1 && k + 1 <= m_wf.hi {
            let k_index = crate::utils::compute_k_index(m_wf.vals.len(), k, m_wf.hi);
						// make safe
						let d_s_k = m_wf.vals[k_index+1];
						Some(d_s_k as isize)
        } else {
            None
        }
		}



		pub fn backtrace_insertion_extend_offset(
				all_wavefronts: &crate::types::WaveFronts,
				score: isize,
				k: i32
		) -> Option<isize> {
				if score < 0 {
						return None;
				}

				// we know it's >= 0
				let score = score as usize;

				// make safe use get to return None
				let i_wf = &all_wavefronts.wavefront_set[score].i;

				if i_wf.lo <= k - 1 && k - 1 <= i_wf.hi {
            let k_index = crate::utils::compute_k_index(i_wf.vals.len(), k, i_wf.hi);
						// make safe
						let i_s_k = i_wf.vals[k_index-1];
						Some(i_s_k as isize)
        } else {
            None
        }
				}

		pub fn backtrace_insertion_open_offset(
				all_wavefronts: &crate::types::WaveFronts,
				score: isize,
				k: i32
		) -> Option<isize> {
				if score < 0 {
						return None;
				}

				// we know it's >= 0
				let score = score as usize;

				// make safe use get to return None
				let m_wf = &all_wavefronts.wavefront_set[score].m;

				if m_wf.lo <= k - 1 && k - 1 <= m_wf.hi {
            let k_index = crate::utils::compute_k_index(m_wf.vals.len(), k, m_wf.hi);
						// make safe
						let i_s_k = m_wf.vals[k_index-1];
						Some(i_s_k as isize + 1)
        } else {
            None
        }
		}

		pub fn backtrace_mismatch_offset(
				all_wavefronts: &crate::types::WaveFronts,
				score: isize,
				k: i32
		) -> Option<isize> {
				if score < 0 {
						return None;
				}

				// we know it's >= 0
				let score = score as usize;

				// make safe use get to return None
				let m_wf = &all_wavefronts.wavefront_set[score].m;

				if m_wf.lo <= k && k <= m_wf.hi {
            let k_index = crate::utils::compute_k_index(m_wf.vals.len(), k, m_wf.hi);
						// make safe
						let m_s_k = m_wf.vals[k_index];
						Some(m_s_k as isize + 1)
        } else {
            None
        }
		}


}
pub mod debug_utils {
		use super::*;

		pub fn viz(v: &Vec<Vec<i32>>) {
				let dim = v[0].len();
				let mut matrix = Array2::<u32>::zeros((dim, dim));

				for i in 0..dim {
						for j in 0..dim {
								matrix[[i,j]] = (v[i][j] + 1) as u32;
						}
				}

				let config = ndarray_to_img::Config {
						verbosity: 0,
						annotate_image: true,
						draw_diagonal: true,
						draw_boundaries: true,
						scaling_factor: 100,
				};

				let scaled_matrix = ndarray_to_img::scale_matrix(&matrix, &config);
				let image_name = format!("m_image.png"); 
				ndarray_to_img::generate_image(&scaled_matrix, &config, &image_name).unwrap();
		}

		pub fn visualize_all(
				all_wavefronts: &types::WaveFronts,
				a_offset: u32,
				each_wf: &Vec<types::WfType>
		) {
				// let dim = all_wavefronts.wavefront_set.len();

				let dim = (a_offset as usize+10, a_offset as usize+10);
				let mut matrix = Array2::<u32>::zeros(dim);
				
				for wf_set in all_wavefronts.wavefront_set.iter() {
						for wf in each_wf {
								let wf_specific = match wf {
										types::WfType::I => &wf_set.i,
										types::WfType::D => &wf_set.d,
										types::WfType::M => &wf_set.m,
								};

								let lo = wf_specific.lo;
								let hi = wf_specific.hi;
								let vals = &wf_specific.vals;
								let len = vals.len();

								for k in lo..=hi {

										let k_index: usize = compute_k_index(len, k, hi);
										let m_s_k: i32 = vals[k_index];

										for offset in  0..=m_s_k {
												let v: usize = compute_v(offset, k);
												let h: usize = compute_h(offset, k);

												// eprintln!("offset: {}\tk: {}\tscore: {}\t({}, {})", m_s_k,  k, score, v, h);

												if v >= dim.0 || h >= dim.0 {
														continue;
												}

												if offset < 0 {
														matrix[[v,h]] = 0;
												} else {
														matrix[[v,h]] = offset as u32 + 1;
												}
										}
								}
						}
				}

				let config = ndarray_to_img::Config {
						verbosity: 0,
						annotate_image: true,
						draw_diagonal: true,
						draw_boundaries: true,
						scaling_factor: 100,
				};

				let scaled_matrix = ndarray_to_img::scale_matrix(&matrix, &config);
				let image_name = "all.png";
				ndarray_to_img::generate_image(&scaled_matrix, &config, &image_name).unwrap();
		}

		pub fn visualize(
				all_wavefronts: &types::WaveFronts,
				a_offset: u32,
				wf_type: types::WfType
		) {
				let a = a_offset as i32;
				let _min_diagonal = -a;

				let dim = all_wavefronts.wavefront_set.len();

				let mut matrix = Array2::<u32>::zeros((dim, dim));

				for (_score, wf) in all_wavefronts.wavefront_set.iter().enumerate() {
						// eprintln!("score: {}   ", score);

						let wf_specific = match wf_type {
								types::WfType::I => &wf.i,
								types::WfType::D => &wf.d,
								types::WfType::M => &wf.m,
						};

						let lo = wf_specific.lo;
						let hi = wf_specific.hi;
						let vals = &wf_specific.vals;
						let len = vals.len();

						for k in lo..=hi {
								let k_index: usize = compute_k_index(len, k, hi);
								let m_s_k: i32 = vals[k_index];

								let v: usize = abs_sub(m_s_k as i32, k) as usize;
								let h: usize = m_s_k as usize;

								// eprintln!("offset: {}\tk: {}\tscore: {}\t({}, {})", m_s_k,  k, score, v, h);

								if v >= dim || h >= dim {
										// eprintln!("\t {:?} ({}, {})", wf_type, v, h);
										continue;
								}

								if m_s_k < 0 {
										matrix[[v,h]] = 0;
								} else {
										matrix[[v,h]] = m_s_k as u32 + 1;
								}
						}
						eprintln!();
				}

				let config = ndarray_to_img::Config {
						verbosity: 0,
						annotate_image: true,
						draw_diagonal: true,
						draw_boundaries: true,
						scaling_factor: 100,
				};

				let scaled_matrix = ndarray_to_img::scale_matrix(&matrix, &config);
				let image_name = format!("{:?}_image.png", wf_type);
				ndarray_to_img::generate_image(&scaled_matrix, &config, &image_name).unwrap();
		}

}

