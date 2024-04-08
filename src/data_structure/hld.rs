use std::ops::{Range, RangeInclusive};

pub struct HLD {
    g: Vec<Vec<usize>>,
    size: Vec<usize>,
    parent: Vec<usize>,
    head: Vec<usize>,
    pub pos: Vec<usize>,
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

    pub fn path(&self, mut u: usize, mut v: usize) -> Vec<RangeInclusive<usize>> {
        let mut res = vec![];
        while self.head[u] != self.head[v] {
            if self.pos[u] > self.pos[v] {
                std::mem::swap(&mut u, &mut v);
            }
            res.push(self.pos[self.head[v]]..=self.pos[v]);
            v = self.parent[self.head[v]];
        }
        if self.pos[u] > self.pos[v] {
            std::mem::swap(&mut u, &mut v);
        }
        res.push(self.pos[u]..=self.pos[v]);
        res
    }

    pub fn subtree(&self, u: usize) -> Range<usize> {
        self.pos[u]..self.pos[u] + self.size[u]
    }
}
