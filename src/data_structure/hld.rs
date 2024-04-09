//! 重軽分解．

use std::ops::{Range, RangeInclusive};

/// 重軽分解．
pub struct HLD {
    g: Vec<Vec<usize>>,
    size: Vec<usize>,
    parent: Vec<usize>,
    head: Vec<usize>,
    pos: Vec<usize>,
}

impl HLD {
    /// 根付き木 $G$，根 $\mathrm{root}$ を受け取り，重軽分解する．
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

    /// $G$ の頂点 $u$ の，$\mathrm{HLD}$ 上の位置を返す．
    pub fn pos(&self, u: usize) -> usize {
        self.pos[u]
    }

    /// $G$ の頂点 $u, v$ の最小共通祖先を返す．
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

    /// $G$ 上の $u$ から $v$ へのパスに対応する，$\mathrm{HLD}$ 上の区間の集合を返す．
    /// 
    /// - `up` は $u$ から $\mathrm{LCA}(u, v)$ へのパスに対応している．
    /// - `down` は $\mathrm{LCA}(u, v)$ から $v$ へのパスに対応している．
    /// 
    /// # Examples
    /// 可換な場合は，そのまま処理する．
    /// ```
    /// // let mut ans = e();
    /// // for range in up {
    /// //     ans = op(ans, seg.prod(range));
    /// // }
    /// // for range in down {
    /// //     ans = op(ans, seg.prod(range));
    /// // }
    /// ```
    /// 
    /// 非可換な場合は，逆順のセグ木を持つなどして処理する．
    /// ```
    /// // let mut ans = e();
    /// // for range in up {
    /// //     ans = op(ans, seg_rev.prod(n - 1 - range.end()..=n - 1 - range.start()));
    /// // }
    /// // for range in down {
    /// //     ans = op(ans, seg.prod(range));
    /// // }
    /// ```
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

    /// $G$ 上の $u$ を根とする部分木に対応する，$\mathrm{HLD}$ 上の区間を返す．
    pub fn subtree(&self, u: usize) -> Range<usize> {
        self.pos[u]..self.pos[u] + self.size[u]
    }
}
