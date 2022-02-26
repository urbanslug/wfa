/*!
Single crate wfa and wflambda.

# Example wflambda

```
use wfa::wflambda;
```

# Example wfa

```
use wfa::wfa;
```


# Visualization

```
use wfa::wflambda;
use ndarray_to_img; // https://github.com/urbanslug/ndarray-to-img

```


Constraints:
- max qlen and tlen is max value of i32 because of
 * hi and lo fields in WaveFrontSet
 * further bound by computing a_k in wf_align
 */

mod core;
mod utils;
mod tests_prelude;

pub mod wfa;
pub mod types;
pub mod wflambda;
