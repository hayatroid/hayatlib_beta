//! 平方剰余

use ac_library::ModInt as M;

/// $X^2 \equiv Y \pmod p$ なる $X$ を返す．
/// そのような $X$ が存在しない場合は `None` を返す．
pub fn sqrt_mod(a: u32, p: u32) -> Option<u32> {
    // Cipolla のアルゴリズム ( https://37zigen.com/cipolla-algorithm/ )
    M::set_modulus(p);

    let a = M::new(a);
    let p = p as u64;
    if p == 2 {
        return Some(a.val());
    }
    if a.val() == 0 {
        return Some(0);
    }
    if a.pow((p - 1) / 2).val() != 1 {
        return None;
    }
    let b: M = (0..).find(|&b| (M::new(b) * b - a).pow((p - 1) / 2).val() != 1).unwrap().into();
    let base = b * b - a;

    let mul = |a: (M, M), b: (M, M)| (a.0 * b.0 + a.1 * b.1 * base, a.0 * b.1 + a.1 * b.0);
    let pow = |mut x: (M, M), mut m: u64| {
        let mut res = (M::new(1), M::new(0));
        while m > 0 {
            if m & 1 == 1 {
                res = mul(res, x);
            }
            x = mul(x, x);
            m >>= 1;
        }
        res
    };

    let t = pow((b, M::new(1)), (p + 1) / 2);
    Some(t.0.val())
}
