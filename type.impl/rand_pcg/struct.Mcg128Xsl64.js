(function() {var type_impls = {
"rand_pcg":[["<details class=\"toggle implementors-toggle\" open><summary><section id=\"impl-Mcg128Xsl64\" class=\"impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#180-222\">source</a><a href=\"#impl-Mcg128Xsl64\" class=\"anchor\">§</a><h3 class=\"code-header\">impl <a class=\"struct\" href=\"rand_pcg/struct.Mcg128Xsl64.html\" title=\"struct rand_pcg::Mcg128Xsl64\">Mcg128Xsl64</a></h3></section></summary><div class=\"impl-items\"><details class=\"toggle method-toggle\" open><summary><section id=\"method.advance\" class=\"method\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#194-211\">source</a><h4 class=\"code-header\">pub fn <a href=\"rand_pcg/struct.Mcg128Xsl64.html#tymethod.advance\" class=\"fn\">advance</a>(&amp;mut self, delta: <a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.u128.html\">u128</a>)</h4></section></summary><div class=\"docblock\"><p>Multi-step advance functions (jump-ahead, jump-back)</p>\n<p>The method used here is based on Brown, “Random Number Generation\nwith Arbitrary Stride,”, Transactions of the American Nuclear\nSociety (Nov. 1994).  The algorithm is very similar to fast\nexponentiation.</p>\n<p>Even though delta is an unsigned integer, we can pass a\nsigned integer to go backwards, it just goes “the long way round”.</p>\n<p>Using this function is equivalent to calling <code>next_64()</code> <code>delta</code>\nnumber of times.</p>\n</div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.new\" class=\"method\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#218-221\">source</a><h4 class=\"code-header\">pub fn <a href=\"rand_pcg/struct.Mcg128Xsl64.html#tymethod.new\" class=\"fn\">new</a>(state: <a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.u128.html\">u128</a>) -&gt; Self</h4></section></summary><div class=\"docblock\"><p>Construct an instance compatible with PCG seed.</p>\n<p>Note that PCG specifies a default value for the parameter:</p>\n<ul>\n<li><code>state = 0xcafef00dd15ea5e5</code></li>\n</ul>\n</div></details></div></details>",0,"rand_pcg::pcg128::Pcg64Mcg"],["<details class=\"toggle implementors-toggle\" open><summary><section id=\"impl-SeedableRng-for-Mcg128Xsl64\" class=\"impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#233-244\">source</a><a href=\"#impl-SeedableRng-for-Mcg128Xsl64\" class=\"anchor\">§</a><h3 class=\"code-header\">impl <a class=\"trait\" href=\"rand_core/trait.SeedableRng.html\" title=\"trait rand_core::SeedableRng\">SeedableRng</a> for <a class=\"struct\" href=\"rand_pcg/struct.Mcg128Xsl64.html\" title=\"struct rand_pcg::Mcg128Xsl64\">Mcg128Xsl64</a></h3></section></summary><div class=\"docblock\"><p>We use a single 126-bit seed to initialise the state and select a stream.\nTwo <code>seed</code> bits (lowest order of last byte) are ignored.</p>\n</div><div class=\"impl-items\"><details class=\"toggle\" open><summary><section id=\"associatedtype.Seed\" class=\"associatedtype trait-impl\"><a href=\"#associatedtype.Seed\" class=\"anchor\">§</a><h4 class=\"code-header\">type <a href=\"rand_core/trait.SeedableRng.html#associatedtype.Seed\" class=\"associatedtype\">Seed</a> = [<a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.u8.html\">u8</a>; <a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.array.html\">16</a>]</h4></section></summary><div class='docblock'>Seed type, which is restricted to types mutably-dereferenceable as <code>u8</code>\narrays (we recommend <code>[u8; N]</code> for some <code>N</code>). <a href=\"rand_core/trait.SeedableRng.html#associatedtype.Seed\">Read more</a></div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.from_seed\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#236-243\">source</a><a href=\"#method.from_seed\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"rand_core/trait.SeedableRng.html#tymethod.from_seed\" class=\"fn\">from_seed</a>(seed: Self::<a class=\"associatedtype\" href=\"rand_core/trait.SeedableRng.html#associatedtype.Seed\" title=\"type rand_core::SeedableRng::Seed\">Seed</a>) -&gt; Self</h4></section></summary><div class='docblock'>Create a new PRNG using the given seed. <a href=\"rand_core/trait.SeedableRng.html#tymethod.from_seed\">Read more</a></div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.seed_from_u64\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/rand_core/lib.rs.html#335\">source</a><a href=\"#method.seed_from_u64\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"rand_core/trait.SeedableRng.html#method.seed_from_u64\" class=\"fn\">seed_from_u64</a>(state: <a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.u64.html\">u64</a>) -&gt; Self</h4></section></summary><div class='docblock'>Create a new PRNG using a <code>u64</code> seed. <a href=\"rand_core/trait.SeedableRng.html#method.seed_from_u64\">Read more</a></div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.from_rng\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/rand_core/lib.rs.html#390\">source</a><a href=\"#method.from_rng\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"rand_core/trait.SeedableRng.html#method.from_rng\" class=\"fn\">from_rng</a>&lt;R&gt;(rng: R) -&gt; <a class=\"enum\" href=\"https://doc.rust-lang.org/1.77.0/core/result/enum.Result.html\" title=\"enum core::result::Result\">Result</a>&lt;Self, <a class=\"struct\" href=\"rand_core/error/struct.Error.html\" title=\"struct rand_core::error::Error\">Error</a>&gt;<div class=\"where\">where\n    R: <a class=\"trait\" href=\"rand_core/trait.RngCore.html\" title=\"trait rand_core::RngCore\">RngCore</a>,</div></h4></section></summary><div class='docblock'>Create a new PRNG seeded from another <code>Rng</code>. <a href=\"rand_core/trait.SeedableRng.html#method.from_rng\">Read more</a></div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.from_entropy\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/rand_core/lib.rs.html#412\">source</a><a href=\"#method.from_entropy\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"rand_core/trait.SeedableRng.html#method.from_entropy\" class=\"fn\">from_entropy</a>() -&gt; Self</h4></section></summary><div class='docblock'>Creates a new instance of the RNG seeded via <a href=\"https://docs.rs/getrandom\"><code>getrandom</code></a>. <a href=\"rand_core/trait.SeedableRng.html#method.from_entropy\">Read more</a></div></details></div></details>","SeedableRng","rand_pcg::pcg128::Pcg64Mcg"],["<details class=\"toggle implementors-toggle\" open><summary><section id=\"impl-PartialEq-for-Mcg128Xsl64\" class=\"impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#171\">source</a><a href=\"#impl-PartialEq-for-Mcg128Xsl64\" class=\"anchor\">§</a><h3 class=\"code-header\">impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/cmp/trait.PartialEq.html\" title=\"trait core::cmp::PartialEq\">PartialEq</a> for <a class=\"struct\" href=\"rand_pcg/struct.Mcg128Xsl64.html\" title=\"struct rand_pcg::Mcg128Xsl64\">Mcg128Xsl64</a></h3></section></summary><div class=\"impl-items\"><details class=\"toggle method-toggle\" open><summary><section id=\"method.eq\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#171\">source</a><a href=\"#method.eq\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"https://doc.rust-lang.org/1.77.0/core/cmp/trait.PartialEq.html#tymethod.eq\" class=\"fn\">eq</a>(&amp;self, other: &amp;<a class=\"struct\" href=\"rand_pcg/struct.Mcg128Xsl64.html\" title=\"struct rand_pcg::Mcg128Xsl64\">Mcg128Xsl64</a>) -&gt; <a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.bool.html\">bool</a></h4></section></summary><div class='docblock'>This method tests for <code>self</code> and <code>other</code> values to be equal, and is used\nby <code>==</code>.</div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.ne\" class=\"method trait-impl\"><span class=\"rightside\"><span class=\"since\" title=\"Stable since Rust version 1.0.0\">1.0.0</span> · <a class=\"src\" href=\"https://doc.rust-lang.org/1.77.0/src/core/cmp.rs.html#242\">source</a></span><a href=\"#method.ne\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"https://doc.rust-lang.org/1.77.0/core/cmp/trait.PartialEq.html#method.ne\" class=\"fn\">ne</a>(&amp;self, other: <a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.reference.html\">&amp;Rhs</a>) -&gt; <a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.bool.html\">bool</a></h4></section></summary><div class='docblock'>This method tests for <code>!=</code>. The default implementation is almost always\nsufficient, and should not be overridden without very good reason.</div></details></div></details>","PartialEq","rand_pcg::pcg128::Pcg64Mcg"],["<details class=\"toggle implementors-toggle\" open><summary><section id=\"impl-Clone-for-Mcg128Xsl64\" class=\"impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#171\">source</a><a href=\"#impl-Clone-for-Mcg128Xsl64\" class=\"anchor\">§</a><h3 class=\"code-header\">impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/clone/trait.Clone.html\" title=\"trait core::clone::Clone\">Clone</a> for <a class=\"struct\" href=\"rand_pcg/struct.Mcg128Xsl64.html\" title=\"struct rand_pcg::Mcg128Xsl64\">Mcg128Xsl64</a></h3></section></summary><div class=\"impl-items\"><details class=\"toggle method-toggle\" open><summary><section id=\"method.clone\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#171\">source</a><a href=\"#method.clone\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"https://doc.rust-lang.org/1.77.0/core/clone/trait.Clone.html#tymethod.clone\" class=\"fn\">clone</a>(&amp;self) -&gt; <a class=\"struct\" href=\"rand_pcg/struct.Mcg128Xsl64.html\" title=\"struct rand_pcg::Mcg128Xsl64\">Mcg128Xsl64</a></h4></section></summary><div class='docblock'>Returns a copy of the value. <a href=\"https://doc.rust-lang.org/1.77.0/core/clone/trait.Clone.html#tymethod.clone\">Read more</a></div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.clone_from\" class=\"method trait-impl\"><span class=\"rightside\"><span class=\"since\" title=\"Stable since Rust version 1.0.0\">1.0.0</span> · <a class=\"src\" href=\"https://doc.rust-lang.org/1.77.0/src/core/clone.rs.html#169\">source</a></span><a href=\"#method.clone_from\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"https://doc.rust-lang.org/1.77.0/core/clone/trait.Clone.html#method.clone_from\" class=\"fn\">clone_from</a>(&amp;mut self, source: <a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.reference.html\">&amp;Self</a>)</h4></section></summary><div class='docblock'>Performs copy-assignment from <code>source</code>. <a href=\"https://doc.rust-lang.org/1.77.0/core/clone/trait.Clone.html#method.clone_from\">Read more</a></div></details></div></details>","Clone","rand_pcg::pcg128::Pcg64Mcg"],["<section id=\"impl-StructuralPartialEq-for-Mcg128Xsl64\" class=\"impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#171\">source</a><a href=\"#impl-StructuralPartialEq-for-Mcg128Xsl64\" class=\"anchor\">§</a><h3 class=\"code-header\">impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/marker/trait.StructuralPartialEq.html\" title=\"trait core::marker::StructuralPartialEq\">StructuralPartialEq</a> for <a class=\"struct\" href=\"rand_pcg/struct.Mcg128Xsl64.html\" title=\"struct rand_pcg::Mcg128Xsl64\">Mcg128Xsl64</a></h3></section>","StructuralPartialEq","rand_pcg::pcg128::Pcg64Mcg"],["<details class=\"toggle implementors-toggle\" open><summary><section id=\"impl-RngCore-for-Mcg128Xsl64\" class=\"impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#246-268\">source</a><a href=\"#impl-RngCore-for-Mcg128Xsl64\" class=\"anchor\">§</a><h3 class=\"code-header\">impl <a class=\"trait\" href=\"rand_core/trait.RngCore.html\" title=\"trait rand_core::RngCore\">RngCore</a> for <a class=\"struct\" href=\"rand_pcg/struct.Mcg128Xsl64.html\" title=\"struct rand_pcg::Mcg128Xsl64\">Mcg128Xsl64</a></h3></section></summary><div class=\"impl-items\"><details class=\"toggle method-toggle\" open><summary><section id=\"method.next_u32\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#248-250\">source</a><a href=\"#method.next_u32\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"rand_core/trait.RngCore.html#tymethod.next_u32\" class=\"fn\">next_u32</a>(&amp;mut self) -&gt; <a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.u32.html\">u32</a></h4></section></summary><div class='docblock'>Return the next random <code>u32</code>. <a href=\"rand_core/trait.RngCore.html#tymethod.next_u32\">Read more</a></div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.next_u64\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#253-256\">source</a><a href=\"#method.next_u64\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"rand_core/trait.RngCore.html#tymethod.next_u64\" class=\"fn\">next_u64</a>(&amp;mut self) -&gt; <a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.u64.html\">u64</a></h4></section></summary><div class='docblock'>Return the next random <code>u64</code>. <a href=\"rand_core/trait.RngCore.html#tymethod.next_u64\">Read more</a></div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.fill_bytes\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#259-261\">source</a><a href=\"#method.fill_bytes\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"rand_core/trait.RngCore.html#tymethod.fill_bytes\" class=\"fn\">fill_bytes</a>(&amp;mut self, dest: &amp;mut [<a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.u8.html\">u8</a>])</h4></section></summary><div class='docblock'>Fill <code>dest</code> with random data. <a href=\"rand_core/trait.RngCore.html#tymethod.fill_bytes\">Read more</a></div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.try_fill_bytes\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#264-267\">source</a><a href=\"#method.try_fill_bytes\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"rand_core/trait.RngCore.html#tymethod.try_fill_bytes\" class=\"fn\">try_fill_bytes</a>(&amp;mut self, dest: &amp;mut [<a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.u8.html\">u8</a>]) -&gt; <a class=\"enum\" href=\"https://doc.rust-lang.org/1.77.0/core/result/enum.Result.html\" title=\"enum core::result::Result\">Result</a>&lt;<a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.unit.html\">()</a>, <a class=\"struct\" href=\"rand_core/error/struct.Error.html\" title=\"struct rand_core::error::Error\">Error</a>&gt;</h4></section></summary><div class='docblock'>Fill <code>dest</code> entirely with random data. <a href=\"rand_core/trait.RngCore.html#tymethod.try_fill_bytes\">Read more</a></div></details></div></details>","RngCore","rand_pcg::pcg128::Pcg64Mcg"],["<section id=\"impl-Eq-for-Mcg128Xsl64\" class=\"impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#171\">source</a><a href=\"#impl-Eq-for-Mcg128Xsl64\" class=\"anchor\">§</a><h3 class=\"code-header\">impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/cmp/trait.Eq.html\" title=\"trait core::cmp::Eq\">Eq</a> for <a class=\"struct\" href=\"rand_pcg/struct.Mcg128Xsl64.html\" title=\"struct rand_pcg::Mcg128Xsl64\">Mcg128Xsl64</a></h3></section>","Eq","rand_pcg::pcg128::Pcg64Mcg"],["<details class=\"toggle implementors-toggle\" open><summary><section id=\"impl-Debug-for-Mcg128Xsl64\" class=\"impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#225-229\">source</a><a href=\"#impl-Debug-for-Mcg128Xsl64\" class=\"anchor\">§</a><h3 class=\"code-header\">impl <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/fmt/trait.Debug.html\" title=\"trait core::fmt::Debug\">Debug</a> for <a class=\"struct\" href=\"rand_pcg/struct.Mcg128Xsl64.html\" title=\"struct rand_pcg::Mcg128Xsl64\">Mcg128Xsl64</a></h3></section></summary><div class=\"impl-items\"><details class=\"toggle method-toggle\" open><summary><section id=\"method.fmt\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/rand_pcg/pcg128.rs.html#226-228\">source</a><a href=\"#method.fmt\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"https://doc.rust-lang.org/1.77.0/core/fmt/trait.Debug.html#tymethod.fmt\" class=\"fn\">fmt</a>(&amp;self, f: &amp;mut <a class=\"struct\" href=\"https://doc.rust-lang.org/1.77.0/core/fmt/struct.Formatter.html\" title=\"struct core::fmt::Formatter\">Formatter</a>&lt;'_&gt;) -&gt; <a class=\"type\" href=\"https://doc.rust-lang.org/1.77.0/core/fmt/type.Result.html\" title=\"type core::fmt::Result\">Result</a></h4></section></summary><div class='docblock'>Formats the value using the given formatter. <a href=\"https://doc.rust-lang.org/1.77.0/core/fmt/trait.Debug.html#tymethod.fmt\">Read more</a></div></details></div></details>","Debug","rand_pcg::pcg128::Pcg64Mcg"]]
};if (window.register_type_impls) {window.register_type_impls(type_impls);} else {window.pending_type_impls = type_impls;}})()