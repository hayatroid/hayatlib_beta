//! 密な形式的冪級数，疎な形式的冪級数，多項式についてのライブラリ．

use ac_library::ModInt998244353 as M;

pub mod fps_impl;

/// 密な形式的冪級数．
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FPS {
    pub coef: Vec<M>
}
