use crate::types;
pub mod backtrace;
pub mod debug;

// How many cells does the wave cross?
pub fn compute_wave_length(lo: i32, hi: i32) -> usize {
    (hi - lo + 1) as usize
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

pub fn end_reached(m_wavefront: Option<&types::WaveFront>, a_k: usize, a_offset: u32) -> bool {
    let m_wavefront = match m_wavefront {
        Some(wf) => wf,
        _ => return false,
    };

    let k_index: usize = compute_k_index(m_wavefront.len(), a_k as i32, m_wavefront.hi);
    let m_s_k = m_wavefront.offsets[k_index]; // m_offset

    m_s_k > 0 && (m_s_k as u32) >= a_offset
}
