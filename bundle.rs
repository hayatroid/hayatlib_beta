// hayatlib {{{
// https://ngtkana.github.io/ac-adapter-rs/hayatlib/index.html

#[allow(dead_code)]
mod hayatlib {
    mod sample {
        pub fn aplusb(a: u64, b: u64) -> u64 {
            a + b
        }
    }
    mod data_structure {
        mod hld {
            use std::ops::{Range, RangeInclusive};
            pub struct HLD {
                g: Vec<Vec<usize>>,
                size: Vec<usize>,
                parent: Vec<usize>,
                head: Vec<usize>,
                pos: Vec<usize>,
            }
            impl HLD {
                pub fn new(g: &Vec<Vec<usize>>, root: usize) -> Self {
                    let n = g.len();
                    let mut res = Self {
                        g: g.clone(),
                        size: vec![1; n],
                        parent: vec![root; n],
                        head: vec![root; n],
                        pos: vec![!0; n],
                    };
                    res.dfs_size(root, root);
                    res.dfs_hld(root, root, &mut 0);
                    res
                }
                fn dfs_size(&mut self, v: usize, p: usize) {
                    let len = self.g[v].len();
                    if len > 0 && self.g[v][0] == p {
                        self.g[v].swap(0, len - 1);
                    }
                    for i in 0..len {
                        let c = self.g[v][i];
                        if c == p {
                            continue;
                        }
                        self.dfs_size(c, v);
                        self.size[v] += self.size[c];
                        if self.size[self.g[v][0]] < self.size[c] {
                            self.g[v].swap(0, i);
                        }
                    }
                }
                fn dfs_hld(&mut self, v: usize, p: usize, cnt: &mut usize) {
                    self.parent[v] = p;
                    self.pos[v] = *cnt;
                    *cnt += 1;
                    for i in 0..self.g[v].len() {
                        let c = self.g[v][i];
                        if c == p {
                            continue;
                        }
                        self.head[c] = if c == self.g[v][0] { self.head[v] } else { c };
                        self.dfs_hld(c, v, cnt);
                    }
                }
                pub fn pos(&self, u: usize) -> usize {
                    self.pos[u]
                }
                pub fn lca(&self, mut u: usize, mut v: usize) -> usize {
                    while self.head[u] != self.head[v] {
                        if self.pos[u] > self.pos[v] {
                            std::mem::swap(&mut u, &mut v);
                        }
                        v = self.parent[self.head[v]];
                    }
                    if self.pos[u] > self.pos[v] {
                        std::mem::swap(&mut u, &mut v);
                    }
                    u
                }
                pub fn path(&self, mut u: usize, mut v: usize) -> (Vec<RangeInclusive<usize>>, Vec<RangeInclusive<usize>>) {
                    let mut up = vec![];
                    let mut down = vec![];
                    while self.head[u] != self.head[v] {
                        if self.pos[u] > self.pos[v] {
                            up.push(self.pos[self.head[u]]..=self.pos[u]);
                            u = self.parent[self.head[u]];
                        } else {
                            down.push(self.pos[self.head[v]]..=self.pos[v]);
                            v = self.parent[self.head[v]];
                        }
                    }
                    if self.pos[u] > self.pos[v] {
                        up.push(self.pos[v]..=self.pos[u]);
                    } else {
                        down.push(self.pos[u]..=self.pos[v]);
                    }
                    down.reverse();
                    (up, down)
                }
                pub fn subtree(&self, u: usize) -> Range<usize> {
                    self.pos[u]..self.pos[u] + self.size[u]
                }
            }
        }
        pub use hld::HLD;
    }
    mod polynomial {
        mod fps {
            use std::{collections::VecDeque, ops::*};
            use ac_library::{convolution, ModInt998244353 as M, RemEuclidU32};
            #[derive(Debug, Clone, PartialEq, Eq)]
            pub struct FPS {
                pub coef: Vec<M>
            }
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
                pub fn new() -> Self {
                    Self::from([0])
                }
                pub fn len(&self) -> usize {
                    self.coef.len()
                }
                pub fn pre(&self, len: usize) -> Self {
                    (0..len).map(|i| if i < self.len() { self[i] } else { M::new(0) }).collect()
                }
                
                pub fn diff(&self) -> Self {
                    (1..self.len()).map(|i| self[i] * i).collect()
                }
                pub fn integral(&self) -> Self {
                    (0..self.len()).map(|i| self[i] / (i + 1)).collect::<Self>() << 1
                }
                pub fn inv(&self, len: usize) -> Self {
                    assert!(self.len() > 0 && self[0].val() != 0);
                    let mut g = Self::from([self[0].inv().val()]);
                    for i in 1..=len.next_power_of_two().trailing_zeros() {
                        // g ← g(2 - fg)
                        g = (g.clone() * (Self::from([2]) - self.pre(1 << i) * g.clone()).pre(1 << i)).pre(1 << i);
                    }
                    g.pre(len)
                }
                pub fn log(&self, len: usize) -> Self {
                    assert!(self.len() > 0 && self[0].val() == 1);
                    (self.diff() * self.inv(len)).pre(len - 1).integral()
                }
                pub fn exp(&self, len: usize) -> Self {
                    assert!(self.len() == 0 || self[0].val() == 0);
                    let mut g = Self::from([1]);
                    for i in 1..=len.next_power_of_two().trailing_zeros() {
                        // g ← g(f + 1 - log(g))
                        g = (g.clone() * (self.pre(1 << i) + Self::from([1]) - g.log(1 << i))).pre(1 << i);
                    }
                    g.pre(len)
                }
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
        }
        pub use fps::FPS;
    }
}
// }}}
