/*!
Combined wfa and wflambda implementations into one crate.

Constraints:
- max qlen and tlen is max value of i32 because of
 * hi and lo fields in WaveFrontSet
 * further bound by computing a_k in wf_align
 */

mod core;
mod utils;

pub mod wfa;
pub mod types;
pub mod wflambda;

mod tests_prelude {
    #[allow(dead_code)]
    pub static TEST_CONFIG: crate::types::Config = crate::types::Config {
        adapt: false,
        verbosity: 0,
        penalties: crate::types::Penalties {
            mismatch: 4,
            matches: 0,
            gap_open: 6,
            gap_extend: 2,
        },
    };
    pub use crate::wflambda::wf_align as wflambda_align;
    pub use crate::wfa::wf_align as wfa_align;
}
