use std::{collections::VecDeque, ops::*};

use ac_library::{convolution, ModInt998244353 as M, RemEuclidU32};

use super::FPS;


impl From<Vec<M>> for FPS {
    fn from(value: Vec<M>) -> Self {
        Self {
            coef: value
        }
    }
}

impl<T: Copy + RemEuclidU32, const N: usize> From<[T; N]> for FPS {
    fn from(value: [T; N]) -> Self {
        Self {
            coef: value.iter().map(|&x| M::new(x)).collect()
        }
    }
}

impl FromIterator<M> for FPS {
    fn from_iter<T: IntoIterator<Item = M>>(iter: T) -> Self {
        Self {
            coef: iter.into_iter().collect()
        }
    }
}

impl Index<usize> for FPS {
    type Output = M;
    fn index(&self, index: usize) -> &Self::Output {
        &self.coef[index]
    }
}

impl IndexMut<usize> for FPS {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.coef[index]
    }
}

impl Add for FPS {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let n = self.len().max(rhs.len());
        (0..n).map(|i| self.coef.get(i).unwrap_or(&M::new(0)) + rhs.coef.get(i).unwrap_or(&M::new(0))).collect()
    }
}

impl Sub for FPS {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let n = self.len().max(rhs.len());
        (0..n).map(|i| self.coef.get(i).unwrap_or(&M::new(0)) - rhs.coef.get(i).unwrap_or(&M::new(0))).collect()
    }
}

impl Mul for FPS {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            coef: convolution(&self.coef, &rhs.coef)
        }
    }
}

impl Mul<M> for FPS {
    type Output = Self;
    fn mul(self, rhs: M) -> Self::Output {
        Self {
            coef: self.coef.iter().map(|&x| x * rhs).collect()
        }
    }
}

impl Div<M> for FPS {
    type Output = Self;
    fn div(self, rhs: M) -> Self::Output {
        let rhs_inv = rhs.inv();
        Self {
            coef: self.coef.iter().map(|&x| x * rhs_inv).collect()
        }
    }
}

impl Shl<usize> for FPS {
    type Output = Self;
    fn shl(self, rhs: usize) -> Self::Output {
        (0..self.len() + rhs).map(|i| if i < rhs { M::new(0) } else { self[i - rhs] }).collect()
    }
}

impl Shr<usize> for FPS {
    type Output = Self;
    fn shr(self, rhs: usize) -> Self::Output {
        (rhs..self.len()).map(|i| self[i]).collect()
    }
}

impl AddAssign for FPS {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}

impl SubAssign for FPS {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.clone() - rhs;
    }
}

impl FPS {

    /// $f(x) = 0$ を返す．
    pub fn new() -> Self {
        Self::from([0])
    }

    /// $f(x)$ の長さを返す．
    pub fn len(&self) -> usize {
        self.coef.len()
    }

    /// $f(x)$ の先頭 $\\mathrm{len}$ 項を返す．
    pub fn pre(&self, len: usize) -> Self {
        (0..len).map(|i| if i < self.len() { self[i] } else { M::new(0) }).collect()
    }
    
    /// $f\'(x)$ を返す．$\\mathrm{len}$ が $1$ 減る．
    /// 
    /// # Examples
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([1, 1, 1, 1, 1]);
    /// let g = FPS::from([1, 2, 3, 4]);
    /// assert_eq!(f.diff(), g);
    /// ```
    pub fn diff(&self) -> Self {
        (1..self.len()).map(|i| self[i] * i).collect()
    }

    /// $\\int_0^x f(x) dx$ を返す．$\\mathrm{len}$ が $1$ 増える．
    /// 
    /// # Examples
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([1, 2, 3, 4]);
    /// let g = FPS::from([0, 1, 1, 1, 1]);
    /// assert_eq!(f.integral(), g);
    /// ```
    pub fn integral(&self) -> Self {
        (0..self.len()).map(|i| self[i] / (i + 1)).collect::<Self>() << 1
    }

    /// $f(x)g(x) = 1$ なる $g(x)$ の先頭 $\\mathrm{len}$ 項を返す．
    /// 
    /// # Panics
    /// $[x^0]f(x) = 0$ のとき panic する．
    /// 
    /// # Examples
    /// verified with [Inv of Formal Power Series](https://judge.yosupo.jp/problem/inv_of_formal_power_series)
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([5, 4, 3, 2, 1]);
    /// let g = FPS::from([598946612, 718735934, 862483121, 635682004, 163871793]);
    /// assert_eq!(f.inv(5), g);
    /// ```
    pub fn inv(&self, len: usize) -> Self {
        assert!(self.len() > 0 && self[0].val() != 0);
        let mut g = Self::from([self[0].inv().val()]);
        for i in 1..=len.next_power_of_two().trailing_zeros() {
            // g ← g(2 - fg)
            g = (g.clone() * (Self::from([2]) - self.pre(1 << i) * g.clone()).pre(1 << i)).pre(1 << i);
        }
        g.pre(len)
    }

    /// $\\int_0^x \\frac{f\'(x)}{f(x)} dx$ の先頭 $\\mathrm{len}$ 項を返す．
    /// 
    /// # Panics
    /// $[x^0]f(x) \\neq 1$ のとき panic する．
    /// 
    /// # Examples
    /// verified with [Log of Formal Power Series](https://judge.yosupo.jp/problem/log_of_formal_power_series)
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([1, 1, 499122179, 166374064, 291154613]);
    /// let g = FPS::from([0, 1, 2, 3, 4]);
    /// assert_eq!(f.log(5), g);
    /// ```
    pub fn log(&self, len: usize) -> Self {
        assert!(self.len() > 0 && self[0].val() == 1);
        (self.diff() * self.inv(len)).pre(len - 1).integral()
    }

    /// $\\log g(x) = f(x)$ なる $g(x)$ の先頭 $\\mathrm{len}$ 項を返す．
    /// 
    /// # Panics
    /// $[x^0]f(x) \\neq 0$ のとき panic する．
    /// 
    /// # Examples
    /// verified with [Exp of Formal Power Series](https://judge.yosupo.jp/problem/exp_of_formal_power_series)
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([0, 1, 2, 3, 4]);
    /// let g = FPS::from([1, 1, 499122179, 166374064, 291154613]);
    /// assert_eq!(f.exp(5), g);
    /// ```
    pub fn exp(&self, len: usize) -> Self {
        assert!(self.len() == 0 || self[0].val() == 0);
        let mut g = Self::from([1]);
        for i in 1..=len.next_power_of_two().trailing_zeros() {
            // g ← g(f + 1 - log(g))
            g = (g.clone() * (self.pre(1 << i) + Self::from([1]) - g.log(1 << i))).pre(1 << i);
        }
        g.pre(len)
    }

    /// $f^m(x)$ の先頭 $\\mathrm{len}$ 項を返す．
    /// 
    /// # Examples
    /// verified with [Pow of Formal Power Series](https://judge.yosupo.jp/problem/pow_of_formal_power_series)
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([0, 0, 9, 12]);
    /// let g = FPS::from([0, 0, 0, 0]);
    /// assert_eq!(f.pow(3, 4), g);
    /// ```
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([1, 1]);
    /// let g = FPS::from([1, 2]);
    /// assert_eq!(f.pow(2, 2), g);
    /// ```
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([0, 0]);
    /// let g = FPS::from([1, 0]);
    /// assert_eq!(f.pow(0, 2), g);
    /// ```
    pub fn pow(&self, m: usize, len: usize) -> Self {
        if m == 0 {
            return Self::from([1]).pre(len);
        }
        if self.coef.iter().all(|&x| x.val() == 0) {
            return Self::from([0]).pre(len);
        }
        // f = a_p x^p f' ([x⁰]f' = 1) となるように f' をとると
        // f^m = (a_p x^p)^m exp(m log(f')) となる
        let (p, &a_p) = self.coef.iter().enumerate().find(|&(_, &x)| x.val() != 0).unwrap();
        if p.saturating_mul(m) >= len {
            return Self::from([0]).pre(len);
        }
        let mut g = (self.clone() >> p) / a_p;
        g = (g.log(len - p * m) * M::new(m)).exp(len - p * m);
        g = (g * a_p.pow(m as u64)) << (p * m);
        g
    }

    /// $g^2(x) = f(x)$ なる $g(x)$ の先頭 $\\mathrm{len}$ 項を返す．
    /// そのような $g(x)$ が存在しない場合は `None` を返す．
    /// 
    /// # Examples
    /// verified with [Sqrt of Formal Power Series](https://judge.yosupo.jp/problem/sqrt_of_formal_power_series)
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([0, 0, 9, 12]);
    /// assert_eq!(f.sqrt(4).unwrap().pow(2, 4), f);
    /// ```
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([0, 0, 10, 12]);
    /// assert_eq!(f.sqrt(4), None);
    /// ```
    pub fn sqrt(&self, len: usize) -> Option<Self> {
        if self.coef.iter().all(|&x| x.val() == 0) {
            return Some(Self::from([0]).pre(len));
        }
        let (p, &a_p) = self.coef.iter().enumerate().find(|&(_, &x)| x.val() != 0).unwrap();
        if p % 2 == 1 {
            return None;
        }
        let f = self.clone() >> p;
        let mut g = Self::from([a_p.sqrt()?.val()]);
        for i in 1..=(len - p / 2).next_power_of_two().trailing_zeros() {
            // g ← (g + f/g) / 2
            g = (g.clone() + (f.pre(1 << i) * g.inv(1 << i)).pre(1 << i)) / M::new(2);
        }
        g = g.pre(len - p / 2) << (p / 2);
        Some(g)
    }

    /// $f \\circ g(x)$ の先頭 $\\mathrm{len}$ 項を返す．
    /// 
    /// # Examples
    /// verified with [Composition of Formal Power Series](https://judge.yosupo.jp/problem/composition_of_formal_power_series)
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([5, 4, 3, 2, 1]);
    /// let g = FPS::from([0, 1, 2, 3, 4]);
    /// let h = FPS::from([5, 4, 11, 26, 59]);
    /// assert_eq!(f.composition(&g, 5), h);
    /// ```
    pub fn composition(&self, rhs: &Self, len: usize) -> Self {
        // k = ⌈ √len ⌉
        let k = (1..).find(|&i| i * i >= len).unwrap();
        // g^0, g^1, ..., g^k を計算
        let mut g = vec![Self::from([1]); k + 1];
        for i in 0..k {
            g[i + 1] = (g[i].clone() * rhs.pre(len)).pre(len);
        }
        let mut res = Self::from([0]);
        for f in self.coef.chunks(k).rev() {
            let mut tmp = Self::from([0]);
            for (i, &f) in f.iter().enumerate() {
                tmp += g[i].clone() * f;
            }
            res = (res * g[k].clone()).pre(len) + tmp;
        }
        res
    }
    
    /// $f \\circ g(x) = x$ なる $g(x)$ の先頭 $\\mathrm{len}$ 項を返す．
    /// # Panics
    /// $[x^0]f(x) \\neq 0$ または $[x^1]f(x) = 0$ のとき panic する．
    /// 
    /// # Examples
    /// verified with [Compositional Inverse of Formal Power Series](https://judge.yosupo.jp/problem/compositional_inverse_of_formal_power_series)
    /// ```
    /// # use hayatlib::polynomial::fps::FPS;
    /// # use ac_library::ModInt998244353 as M;
    /// let f = FPS::from([0, 1, 2, 3, 4]);
    /// let g = FPS::from([0, 1, 998244351, 5, 998244339]);
    /// assert_eq!(f.compositional_inv(5), g);
    /// ```
    pub fn compositional_inv(&self, len: usize) -> Self {
        assert!(self.len() > 1 && self[0].val() == 0 && self[1].val() != 0);
        let mut g = Self::from([0, self[1].inv().val()]);
        for i in 2..=len.next_power_of_two().trailing_zeros() {
            let fg = self.pre(1 << i).composition(&g, 1 << i);
            let dfg = self.pre(1 << i).diff().composition(&g, 1 << i);
            // g ← g - (fg - x) / (f'g)
            g -= ((fg - Self::from([0, 1])) * dfg.inv(1 << i)).pre(1 << i);
        }
        g.pre(len)
    }

    pub fn all_prod(v: &Vec<Self>) -> Self {
        // FFT マージテク ( https://trap.jp/post/1657/ )
        let mut deque = v.iter().cloned().collect::<VecDeque<Self>>();
        deque.push_back(Self::from([1]));
        while deque.len() > 1 {
            let a = deque.pop_front().unwrap();
            let b = deque.pop_front().unwrap();
            deque.push_back(a * b);
        }
        deque.pop_front().unwrap()
    }

    pub fn shrinked(&self) -> Self {
        let mut res = self.clone();
        while res.len() > 0 && res.coef.last().unwrap().val() == 0 {
            res.coef.pop();
        }
        res
    }

    pub fn rev(&self) -> Self {
        self.coef.iter().cloned().rev().collect()
    }

    pub fn rem_euclid(&self, rhs: &Self) -> Self {
        // Euclidean division ( https://cp-algorithms.com/algebra/polynomial.html )
        let a = self.shrinked();
        let b = rhs.shrinked();
        let n = a.len();
        let m = b.len();
        if n < m {
            return a;
        }
        let len = n - m + 1;
        let c = (a.rev().pre(len) * b.rev().inv(len)).pre(len).rev();
        (a - b * c).pre(m - 1)
    }

    pub fn multipoint_evaluation(&self, x: &Vec<M>) -> Vec<M> {
        let n = x.len().next_power_of_two();
        let mut tree = vec![Self::from([1]); n * 2];
        for (i, &x) in x.iter().enumerate() {
            tree[n + i] = Self::from([(-x).val(), 1]);
        }
        for i in (1..n).rev() {
            tree[i] = tree[i * 2].clone() * tree[i * 2 + 1].clone();
        }
        let mut res = vec![Self::from([1]); n * 2];
        res[1] = self.rem_euclid(&tree[1]);
        for i in 2..n * 2 {
            res[i] = res[i / 2].rem_euclid(&tree[i]);
        }
        (0..x.len()).map(|i| res[n + i].coef[0]).collect()
    }

    pub fn multipoint_evaluation_on_geometric_sequence(&self, a: M, r: M, m: usize) -> Vec<M> {
        // Chirp-z Transform ( https://cp-algorithms.com/algebra/polynomial.html )
        if r.val() == 0 {
            let mut res = vec![self[0]; m];
            res[0] = self.coef.iter().rev().fold(M::new(0), |s, &f| s * a + f);
            return res;
        }
        let mut f = self.shrinked();
        let n = f.len();
        if n == 0 {
            return vec![self[0]; m];
        }
        // a = 1 に帰着
        let mut pow = M::new(1);
        for f in f.coef.iter_mut() {
            *f *= pow;
            pow *= a;
        }
        let r_inv = r.inv();
        let mut pow = (M::new(1), M::new(1));
        for f in f.coef.iter_mut() {
            *f *= pow.0;
            pow = (pow.0 * pow.1, pow.1 * r_inv);
        }
        f = f.rev();
        let mut g = Self::from(vec![M::new(1); n + m - 1]);
        let mut pow = (M::new(1), M::new(1));
        for g in g.coef.iter_mut() {
            *g *= pow.0;
            pow = (pow.0 * pow.1, pow.1 * r);
        }
        let mut h = (f * g).coef[n - 1..n + m - 1].iter().cloned().collect::<Self>();
        let mut pow = (M::new(1), M::new(1));
        for h in h.coef.iter_mut() {
            *h *= pow.0;
            pow = (pow.0 * pow.1, pow.1 * r_inv);
        }
        h.coef
    }

    pub fn polynomial_interpolation(x: &Vec<M>, y: &Vec<M>) -> Self {
        assert!(x.len() == y.len());
        // 多項式補間 ( https://37zigen.com/lagrange-interpolation/ )
        let n = x.len().next_power_of_two();
        let mut tree = vec![Self::from([1]); n * 2];
        for (i, &x) in x.iter().enumerate() {
            tree[n + i] = Self::from([(-x).val(), 1]);
        }
        for i in (1..n).rev() {
            tree[i] = tree[i * 2].clone() * tree[i * 2 + 1].clone();
        }
        let mut res = vec![Self::from([0]); n * 2];
        let g_diff = tree[1].diff().multipoint_evaluation(x);
        for (i, &y) in y.iter().enumerate() {
            res[n + i] = Self::from([(y / g_diff[i]).val()]);
        }
        for i in (1..n).rev() {
            res[i] = res[i * 2].clone() * tree[i * 2 + 1].clone() + res[i * 2 + 1].clone() * tree[i * 2].clone();
        }
        res[1].pre(x.len())
    }

    pub fn polynomial_interpolation_on_geometric_sequence(a: M, r: M, y: &Vec<M>) -> Self {
        // 多項式補間 ( https://37zigen.com/lagrange-interpolation/ )
        let m = y.len();
        if m <= 1 {
            return y.clone().into();
        }

        // g = Π (x - r^i) の計算は q-二項係数を使って Θ(n) で出来るらしい
        // 要改善
        let mut v = vec![];
        let mut pow = M::new(1);
        for _ in 0..m {
            v.push(Self::from([(-pow).val(), 1]));
            pow *= r;
        }
        let g = Self::all_prod(&v);

        let g_diff = g.diff().multipoint_evaluation_on_geometric_sequence(M::new(1), r, m);
        let mut h = (0..m).map(|i| -y[i] / g_diff[i]).collect::<Self>();
        h = h.multipoint_evaluation_on_geometric_sequence(r.inv(), r.inv(), m).into();
        let mut res = (g * h).pre(m);

        // r^i で補完してたので，ar^i で補完するように
        let a_inv = a.inv();
        let mut pow = M::new(1);
        for res in res.coef.iter_mut() {
            *res = *res * pow;
            pow *= a_inv;
        }
        res
    }

    pub fn polynomial_taylor_shift(&self, c: M) -> Self {
        // https://twitter.com/risujiroh/status/1215710785000751104
        let mut f = self.shrinked();
        let n = f.len();
        let mut fact = M::new(1);
        for i in 1..n {
            fact *= i;
            f[i] *= fact;
        }
        fact = fact.inv();
        // exp(cx) は線形時間で計算できるらしいですよ
        let mut g = Self::from(vec![M::new(1); n]);
        for i in 1..n {
            g[i] = g[i - 1] * c / i;
        }
        f = (f.rev() * g).pre(n).rev();
        for i in (1..n).rev() {
            f[i] *= fact;
            fact *= i;
        }
        f.pre(self.len())
    }

    pub fn shift_of_sampling_points_of_polynomial(y: &Vec<M>, c: M, m: usize) -> Vec<M> {
        // https://mathlog.info/articles/3351
        let n = y.len();

        // 前計算
        let k = n.max(m);
        let mut fact = vec![M::new(1); k + 1];
        for i in 1..=k {
            fact[i] = fact[i - 1] * i;
        }
        let mut fact_inv = vec![M::new(1); k + 1];
        fact_inv[k] = fact[k].inv();
        for i in (1..=k).rev() {
            fact_inv[i - 1] = fact_inv[i] * i;
        }
        let mut fact_c = vec![M::new(1); k + 1];
        for i in 0..k {
            fact_c[i + 1] = fact_c[i] * (c - i);
        }

        // f の下降冪表現を得る
        let x = (0..n).map(|i| y[i] * fact_inv[i]).collect::<Self>();
        let y = (0..n).map(|i| fact_inv[i] * if i % 2 == 0 { 1 } else { -1 }).collect::<Self>();
        let f = (x * y).pre(n);

        // f のシフトの下降冪表現を得る
        let x = (0..n).map(|i| f[i] * fact[i]).collect::<Self>();
        let y = (0..n).map(|i| fact_c[i] * fact_inv[i]).collect::<Self>().rev();
        let z = x * y;
        let f = (0..n).map(|i| z[n - 1 + i] * fact_inv[i]).collect::<Self>();

        // f(c + i) の計算
        let x = f.clone();
        let y = (0..m).map(|i| fact_inv[i]).collect::<Self>();
        let z = x * y;
        (0..m).map(|i| z[i] * fact[i]).collect()
    }
}

trait Mod998Traits {
    fn sqrt(&self) -> Option<M>;
}

impl Mod998Traits for M {
    // Cipolla のアルゴリズム ( https://37zigen.com/cipolla-algorithm/ )
    fn sqrt(&self) -> Option<M> {
        if self.val() == 0 {
            return Some(M::new(0));
        }
        let p = M::modulus() as u64;
        if self.pow((p - 1) / 2).val() != 1 {
            return None;
        }
        let b: M = (0..).find(|&i| (M::new(i) * i - self).pow((p - 1) / 2).val() != 1).unwrap().into();
        let base = b * b - self;
        let mul = |a: (M, M), b: (M, M)| -> (M, M) {
            (
                a.0 * b.0 + a.1 * b.1 * base,
                a.0 * b.1 + a.1 * b.0,
            )
        };
        let pow = |mut a: (M, M), mut m: u64| -> (M, M) {
            let mut res = (M::new(1), M::new(0));
            while m > 0 {
                if m & 1 == 1 {
                    res = mul(res, a);
                }
                a = mul(a, a);
                m >>= 1;
            }
            res
        };
        Some(pow((b, M::new(1)), (p + 1) / 2).0)
    }
}
