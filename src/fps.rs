//! 形式的冪級数．

use std::ops::*;

use ac_library::{convolution, ModInt998244353 as M, RemEuclidU32};

fn cipolla(a: M) -> M {
    assert!(a.pow((M::modulus() as u64 - 1) / 2).val() == 1);
    if a.val() == 0 {
        return a;
    }
    let mut b = M::new(0);
    while (b * b - a).pow((M::modulus() as u64 - 1) / 2).val() == 1 {
        b += 1;
    }
    let mut n = (M::modulus() + 1) / 2;
    let base = b * b - a;
    let mul = |p: (M, M), q: (M, M)| -> (M, M) {
        (p.0 * q.0 + base * p.1 * q.1, p.0 * q.1 + p.1 * q.0)
    };
    let mut res = (M::new(1), M::new(0));
    let mut p = (b, M::new(1));
    while n > 0 {
        if n & 1 == 1 {
            res = mul(res, p);
        }
        p = mul(p, p);
        n >>= 1;
    }
    res.0
}

/// 形式的冪級数．
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FPS {
    pub coef: Vec<M>,
}

impl FPS {
    /// $f(x) = \\sum \\mathrm{coef}_i x^i$ を返す．
    pub fn new(coef: Vec<M>) -> Self {
        Self { coef }
    }

    /// $f(x) = c$ を返す．
    pub fn constant<T: RemEuclidU32>(c: T) -> Self {
        Self::new(vec![M::new(c)])
    }

    /// $f(x)$ の $\\mathrm{deg}$ を返す．
    pub fn len(&self) -> usize {
        self.coef.len()
    }

    /// $f(x)$ の先頭 $\\mathrm{deg}$ 項を返す．
    pub fn pre(&self, deg: usize) -> Self {
        let mut res = self.clone();
        res.coef.resize(deg, M::new(0));
        res
    }

    /// $f\'(x)$ を返す．$\\mathrm{deg}$ が $1$ 減る．
    /// 
    /// # Examples
    /// ```
    /// # use hayatlib_beta::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([1, 1, 1, 1, 1]);
    /// let g = FPS::from([1, 2, 3, 4]);
    /// assert_eq!(f.diff(), g);
    /// ```
    pub fn diff(&self) -> Self {
        Self::new((1..self.len()).map(|i| self[i] * i).collect())
    }

    /// $\\int_0^x f(t) \\, \\mathrm{d}t$ を返す．$\\mathrm{deg}$ が $1$ 増える．
    /// 
    /// # Examples
    /// ```
    /// # use hayatlib_beta::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([1, 2, 3, 4]);
    /// let g = FPS::from([0, 1, 1, 1, 1]);
    /// assert_eq!(f.integral(), g);
    /// ```
    pub fn integral(&self) -> Self {
        Self::new((0..=self.len()).map(|i| if i > 0 { self[i - 1] * M::new(i).inv() } else { M::new(0) }).collect())
    }

    /// $f(x)g(x) = 1$ なる $g(x)$ の先頭 $\\mathrm{deg}$ 項を返す．
    /// 
    /// # Examples
    /// 参考：[Inv of Formal Power Series](https://judge.yosupo.jp/problem/inv_of_formal_power_series)
    /// ```
    /// # use hayatlib_beta::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([5, 4, 3, 2, 1]);
    /// let g = FPS::from([598946612, 718735934, 862483121, 635682004, 163871793]);
    /// assert_eq!(f.inv(5), g);
    /// ```
    pub fn inv(&self, deg: usize) -> Self {
        assert!(self.coef.get(0).unwrap_or(&M::new(0)).val() > 0);
        let mut res = Self::constant(self[0].inv().val());
        for i in 1..=deg.next_power_of_two().trailing_zeros() {
            // g ← g(2 - fg)
            res = (res.clone() * (Self::constant(2) - self.pre(1 << i) * res.clone())).pre(1 << i);
        }
        res.pre(deg)
    }

    /// $\\int_0^x \\frac{f\'(t)}{f(t)} \\mathrm{d}t$ の先頭 $\\mathrm{deg}$ 項を返す．
    /// 
    /// # Examples
    /// 参考：[Log of Formal Power Series](https://judge.yosupo.jp/problem/log_of_formal_power_series)
    /// ```
    /// # use hayatlib_beta::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([1, 1, 499122179, 166374064, 291154613]);
    /// let g = FPS::from([0, 1, 2, 3, 4]);
    /// assert_eq!(f.log(5), g);
    /// ```
    pub fn log(&self, deg: usize) -> Self {
        assert!(self.coef.get(0).unwrap_or(&M::new(0)).val() == 1);
        (self.diff() * self.inv(deg)).pre(deg - 1).integral()
    }

    /// $\\log g(x) = f(x)$ なる $g(x)$ の先頭 $\\mathrm{deg}$ 項を返す．
    /// 
    /// # Examples
    /// 参考：[Exp of Formal Power Series](https://judge.yosupo.jp/problem/exp_of_formal_power_series)
    /// ```
    /// # use hayatlib_beta::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([0, 1, 2, 3, 4]);
    /// let g = FPS::from([1, 1, 499122179, 166374064, 291154613]);
    /// assert_eq!(f.exp(5), g);
    /// ```
    pub fn exp(&self, deg: usize) -> Self {
        assert!(self.coef.get(0).unwrap_or(&M::new(0)).val() == 0);
        let mut res = Self::constant(1);
        for i in 1..=deg.next_power_of_two().trailing_zeros() {
            // g ← g(1 - log(g) + f)
            res = (res.clone() * (Self::constant(1) - res.log(1 << i) + self.pre(1 << i))).pre(1 << i)
        }
        res.pre(deg)
    }

    /// $f(x)^m$ の先頭 $\\mathrm{deg}$ 項を返す．
    /// 
    /// # Examples
    /// 参考：[Pow of Formal Power Series](https://judge.yosupo.jp/problem/pow_of_formal_power_series)
    /// ```
    /// # use hayatlib_beta::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([0, 0, 9, 12]);
    /// let g = FPS::from([0, 0, 0, 0]);
    /// assert_eq!(f.pow(3, 4), g);
    /// ```
    /// 
    /// ```
    /// # use hayatlib_beta::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([1, 1]);
    /// let g = FPS::from([1, 2]);
    /// assert_eq!(f.pow(2, 2), g);
    /// ```
    /// 
    /// ```
    /// # use hayatlib_beta::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([0, 0]);
    /// let g = FPS::from([1, 0]);
    /// assert_eq!(f.pow(0, 2), g);
    /// ```
    pub fn pow(&self, m: usize, deg: usize) -> Self {
        if m == 0 {
            return FPS::constant(1).pre(deg);
        }
        if self.coef.iter().all(|x| x.val() == 0) {
            return FPS::constant(0).pre(deg);
        }
        let (p, &a_p) = self.coef.iter().enumerate().find(|(_, x)| x.val() > 0).unwrap();
        if p.saturating_mul(m) >= deg {
            return FPS::constant(0).pre(deg);
        }
        let mut res = (self.clone() >> p) * a_p.inv();
        // g^m = e^{m log(g)}
        res = (res.log(deg - p * m) * M::new(m)).exp(deg - p * m);
        res = (res << p * m) * a_p.pow(m as u64);
        res
    }

    /// $g(x)^2 = f(x)$ なる $g(x)$ の先頭 $\\mathrm{deg}$ 項を返す．
    /// 
    /// # Examples
    /// 参考：[Sqrt of Formal Power Series](https://judge.yosupo.jp/problem/sqrt_of_formal_power_series)
    /// ```
    /// # use hayatlib_beta::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([0, 0, 9, 12]);
    /// assert_eq!(f.sqrt(4).unwrap().pow(2, 4), f);
    /// ```
    /// 
    /// ```
    /// # use hayatlib_beta::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([0, 0, 10, 12]);
    /// assert_eq!(f.sqrt(4), None);
    /// ```
    pub fn sqrt(&self, deg: usize) -> Option<Self> {
        if self.coef.iter().all(|x| x.val() == 0) {
            return Some(FPS::constant(0).pre(deg));
        }
        let (p, &a_p) = self.coef.iter().enumerate().find(|(_, x)| x.val() > 0).unwrap();
        if p % 2 == 1 || a_p.pow((M::modulus() as u64 - 1) / 2).val() != 1 {
            return None;
        }
        let f = (self.clone() >> p) * a_p.inv();
        let mut res = Self::constant(1);
        for i in 1..=(deg - p / 2).next_power_of_two().trailing_zeros() {
            // g ← (g + f/g) / 2
            res = (res.clone() + f.pre(1 << i) * res.inv(1 << i)).pre(1 << i) * M::new(2).inv();
        }
        res = (res.pre(deg - p / 2) << p / 2) * cipolla(a_p);
        Some(res)
    }
}

impl<T: Copy + RemEuclidU32, const N: usize> From<[T; N]> for FPS {
    fn from(value: [T; N]) -> Self {
        FPS::new(value.iter().map(|&x| M::new(x)).collect())
    }
}
impl Index<usize> for FPS {
    type Output = M;
    fn index(&self, index: usize) -> &Self::Output {
        &self.coef[index]
    }
}
impl Add for FPS {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self.coef.resize(self.len().max(rhs.len()), M::new(0));
        self.coef.iter_mut().zip(rhs.coef.iter()).for_each(|(a, b)| *a += *b);
        self
    }
}
impl Sub for FPS {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self.coef.resize(self.len().max(rhs.len()), M::new(0));
        self.coef.iter_mut().zip(rhs.coef.iter()).for_each(|(a, b)| *a -= b);
        self
    }
}
impl Mul for FPS {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        FPS::new(convolution(&self.coef, &rhs.coef))
    }
}
impl Mul<M> for FPS {
    type Output = Self;
    fn mul(mut self, rhs: M) -> Self::Output {
        self.coef.iter_mut().for_each(|a| *a *= rhs);
        self
    }
}
impl Shl<usize> for FPS {
    type Output = FPS;
    fn shl(self, rhs: usize) -> Self::Output {
        FPS::new(vec![M::new(0); rhs].iter().chain(self.coef.iter()).cloned().collect())
    }
}
impl Shr<usize> for FPS {
    type Output = FPS;
    fn shr(self, rhs: usize) -> Self::Output {
        FPS::new(self.coef[rhs..].iter().cloned().collect())
    }
}