use num;
use ndarray::{Array2, Array};
use crate::types;

use ndarray_to_img;

// How many cells does the wave cross?
pub fn compute_wave_length(lo: i32, hi: i32) -> usize {
		(hi-lo+1) as usize
}

pub fn to_usize_or_zero<T: num::cast::ToPrimitive>(n: T) -> usize {
		n.to_usize().unwrap_or(0)
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

pub fn end_reached(m_wavefront: Option<&types::WaveFront>, a_k: usize, a_offset: u32) -> bool {
    let m_wavefront = match m_wavefront {
        Some(wf) => wf,
        _ => { return false }
    };

    let k_index: usize = compute_k_index(m_wavefront.len(), a_k as i32, m_wavefront.hi);
    let m_s_k = m_wavefront.offsets[k_index]; // m_offset

    m_s_k >= 0 && (m_s_k as u32) >= a_offset

}


pub mod debug_utils {
    use super::*;

    pub fn visualize_all(
        all_wavefronts: &types::WaveFronts,
        a_offset: u32,
        each_wf: &Vec<types::WfType>
    ) {
        // let dim = all_wavefronts.wavefront_set.len();

        let dim = (a_offset as usize+10, a_offset as usize+10);
        let x = ndarray::Dim(dim);
        let mut matrix: Array2<Option<i32>> = Array::from_elem(x, None);

        for wf_set in all_wavefronts.wavefront_set.iter() {

            let wf_set = match wf_set {
                Some(w) => w,
                None => { continue; },
            };

            for wf in each_wf {
                let wf_specific = match wf {
                    types::WfType::I => &wf_set.i,
                    types::WfType::D => &wf_set.d,
                    types::WfType::M => &wf_set.m,
                };

                let lo = wf_specific.lo;
                let hi = wf_specific.hi;
                let vals = &wf_specific.offsets;
                let len = vals.len();

                for k in lo..=hi {

                    let k_index: usize = compute_k_index(len, k, hi);
                    let m_s_k: i32 = vals[k_index];

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
            scaling_factor: 100,
        };

        let scaled_matrix = ndarray_to_img::scale_matrix(&matrix, &config);
        let image_name = "all.png";
        ndarray_to_img::generate_image(&scaled_matrix, &config, &image_name).unwrap();
    }
}

