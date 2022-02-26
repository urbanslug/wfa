use std::iter;
use std::str;

use num;
use crate::types;
pub mod backtrace;
pub mod macros;

const ASCII_ZERO: u8 = 48;

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

// TODO: make it a macro?
pub fn repeat_char(c: char, count: u32) -> std::iter::Take<std::iter::Repeat<char>> {
    iter::repeat(c).take(count as usize)
}

pub fn vec_u8_to_str_unsafe(v: &Vec<u8>) -> &str {
    str::from_utf8(v).unwrap()
}

// rename to ASCII
pub fn unsigned_literal_to_u8<T: num::Unsigned + num::ToPrimitive>(v: T) -> u8 {
    v.to_u8().unwrap() + ASCII_ZERO
}

// is this a good idea?
pub fn unsigned_num_to_ASCII<T: num::Unsigned + std::fmt::Display>(n: T) -> Vec<u8> {
    n.to_string().as_bytes().to_vec()
}



// Compare current to next and accumulate counts
// O(n)
// RLE for short
pub fn run_length_encode(cigar: &[u8]) -> Vec<u8> {
    match cigar {
        [] => { Vec::<u8>::new() },
        [c] => { Vec::from([ unsigned_literal_to_u8(1 as u8), *c]) },
        [start_char, the_rest@..] => {
            // more than one value
            let mut rle = Vec::<u8>::new();
            let mut current: u8 = *start_char;
            let mut count: u32 = 1;

            let mut update_rle = |count: u32, c: u8| {
                rle.extend_from_slice(&unsigned_num_to_ASCII(count));
                rle.extend_from_slice(&[c]);
            };

            for c in the_rest {
                if current == *c {
                    count += 1;
                } else {
                    update_rle(count, current);

                    // reset
                    current = *c;
                    count = 1;
                }
            }

            // last run
            update_rle(count, current);

            rle
        }
    }
}

// O(n)
pub fn print_aln(cigar: &[u8], t: &[u8], q: &[u8]) {
    let qlen = q.len();
    let tlen = t.len();

    let longer = std::cmp::max(tlen, qlen);
    let mut marker = Vec::<u8>::with_capacity(longer);
    let mut query = Vec::<u8>::with_capacity(qlen);
    let mut text = Vec::<u8>::with_capacity(tlen);

    let space = b' ';
    let dash = b'-';
    let vertical_bar = b'|';

    let mut q_iter = q.iter();
    let mut t_iter = t.iter();

    for c in cigar {
        match c {
            b'M' => {
                query.push(*q_iter.next().unwrap());
                marker.push(vertical_bar);
                text.push(*t_iter.next().unwrap());
            }

            b'X' => {
                query.push(*q_iter.next().unwrap());
                marker.push(space);
                text.push(*t_iter.next().unwrap());
            }

            b'I' => {
                query.push(dash);
                marker.push(space);
                text.push(*t_iter.next().unwrap());
            }

            b'D' => {
                query.push(*q_iter.next().unwrap());
                marker.push(space);
                text.push(dash);
            }

            _ => {
                panic!("[utils::backtrace_utils::print_aln] found char not M, I, X or D")
            }
        }
    }

    eprintln!();
    eprintln!("{}", vec_u8_to_str_unsafe(&query));
    eprintln!("{}", vec_u8_to_str_unsafe(&marker));
    eprintln!("{}", vec_u8_to_str_unsafe(&text));
}


#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::{assert_eq};

    #[test]
    fn test_u32_literal_to_u8() {
        assert_eq!(unsigned_literal_to_u8(0 as u32), 48);
        assert_eq!(unsigned_literal_to_u8(9 as u32), 57);
    }

    #[test]
    fn test_run_length_encode() {
        let cigar = repeat_char('M', 11).map(|x| x as u8 ).collect::<Vec<u8>>();
        assert_eq!(run_length_encode(&cigar), "11M".as_bytes());

        let cigar: Vec<u8> = Vec::from([
            repeat_char('D', 1).map(|x| x as u8 ).collect::<Vec<u8>>(),
            repeat_char('M', 21).map(|x| x as u8 ).collect::<Vec<u8>>(),
            repeat_char('I', 2).map(|x| x as u8 ).collect::<Vec<u8>>(),
            repeat_char('M', 15).map(|x| x as u8 ).collect::<Vec<u8>>(),
            repeat_char('I', 22).map(|x| x as u8 ).collect::<Vec<u8>>(),
            repeat_char('D', 16).map(|x| x as u8 ).collect::<Vec<u8>>(),
        ]).into_iter().flatten().collect();
        assert_eq!(run_length_encode(&cigar), "1D21M2I15M22I16D".as_bytes());
    }
}
