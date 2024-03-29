(function() {var type_impls = {
"ndarray":[["<details class=\"toggle implementors-toggle\" open><summary><section id=\"impl-ArrayBase%3COwnedRepr%3CA%3E,+Dim%3C%5Busize;+0%5D%3E%3E\" class=\"impl\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#20-52\">source</a><a href=\"#impl-ArrayBase%3COwnedRepr%3CA%3E,+Dim%3C%5Busize;+0%5D%3E%3E\" class=\"anchor\">§</a><h3 class=\"code-header\">impl&lt;A&gt; <a class=\"type\" href=\"ndarray/type.Array.html\" title=\"type ndarray::Array\">Array</a>&lt;A, <a class=\"type\" href=\"ndarray/type.Ix0.html\" title=\"type ndarray::Ix0\">Ix0</a>&gt;</h3></section></summary><div class=\"docblock\"><p>Methods specific to <code>Array0</code>.</p>\n<p><em><strong>See also all methods for <a href=\"ndarray/struct.ArrayBase.html\" title=\"struct ndarray::ArrayBase\"><code>ArrayBase</code></a></strong></em></p>\n</div><div class=\"impl-items\"><details class=\"toggle method-toggle\" open><summary><section id=\"method.into_scalar\" class=\"method\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#34-51\">source</a><h4 class=\"code-header\">pub fn <a href=\"ndarray/type.Array.html#tymethod.into_scalar\" class=\"fn\">into_scalar</a>(self) -&gt; A</h4></section></summary><div class=\"docblock\"><p>Returns the single element in the array without cloning it.</p>\n\n<div class=\"example-wrap\"><pre class=\"rust rust-example-rendered\"><code><span class=\"kw\">use </span>ndarray::{arr0, Array0};\n\n<span class=\"comment\">// `Foo` doesn't implement `Clone`.\n</span><span class=\"attr\">#[derive(Debug, Eq, PartialEq)]\n</span><span class=\"kw\">struct </span>Foo;\n\n<span class=\"kw\">let </span>array: Array0&lt;Foo&gt; = arr0(Foo);\n<span class=\"kw\">let </span>scalar: Foo = array.into_scalar();\n<span class=\"macro\">assert_eq!</span>(scalar, Foo);</code></pre></div>\n</div></details></div></details>",0,"ndarray::aliases::Array0"],["<details class=\"toggle implementors-toggle\" open><summary><section id=\"impl-ArrayBase%3COwnedRepr%3CA%3E,+D%3E\" class=\"impl\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#57-69\">source</a><a href=\"#impl-ArrayBase%3COwnedRepr%3CA%3E,+D%3E\" class=\"anchor\">§</a><h3 class=\"code-header\">impl&lt;A, D&gt; <a class=\"type\" href=\"ndarray/type.Array.html\" title=\"type ndarray::Array\">Array</a>&lt;A, D&gt;<div class=\"where\">where\n    D: <a class=\"trait\" href=\"ndarray/trait.Dimension.html\" title=\"trait ndarray::Dimension\">Dimension</a>,</div></h3></section></summary><div class=\"docblock\"><p>Methods specific to <code>Array</code>.</p>\n<p><em><strong>See also all methods for <a href=\"ndarray/struct.ArrayBase.html\" title=\"struct ndarray::ArrayBase\"><code>ArrayBase</code></a></strong></em></p>\n</div><div class=\"impl-items\"><details class=\"toggle method-toggle\" open><summary><section id=\"method.into_raw_vec\" class=\"method\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#66-68\">source</a><h4 class=\"code-header\">pub fn <a href=\"ndarray/type.Array.html#tymethod.into_raw_vec\" class=\"fn\">into_raw_vec</a>(self) -&gt; <a class=\"struct\" href=\"https://doc.rust-lang.org/1.77.0/alloc/vec/struct.Vec.html\" title=\"struct alloc::vec::Vec\">Vec</a>&lt;A&gt;</h4></section></summary><div class=\"docblock\"><p>Return a vector of the elements in the array, in the way they are\nstored internally.</p>\n<p>If the array is in standard memory layout, the logical element order\nof the array (<code>.iter()</code> order) and of the returned vector will be the same.</p>\n</div></details></div></details>",0,"ndarray::aliases::Array0","ndarray::aliases::Array1","ndarray::aliases::Array2","ndarray::aliases::Array3","ndarray::aliases::Array4","ndarray::aliases::Array5","ndarray::aliases::Array6","ndarray::aliases::ArrayD"],["<details class=\"toggle implementors-toggle\" open><summary><section id=\"impl-ArrayBase%3COwnedRepr%3CA%3E,+Dim%3C%5Busize;+2%5D%3E%3E\" class=\"impl\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#74-166\">source</a><a href=\"#impl-ArrayBase%3COwnedRepr%3CA%3E,+Dim%3C%5Busize;+2%5D%3E%3E\" class=\"anchor\">§</a><h3 class=\"code-header\">impl&lt;A&gt; <a class=\"type\" href=\"ndarray/type.Array.html\" title=\"type ndarray::Array\">Array</a>&lt;A, <a class=\"type\" href=\"ndarray/type.Ix2.html\" title=\"type ndarray::Ix2\">Ix2</a>&gt;</h3></section></summary><div class=\"docblock\"><p>Methods specific to <code>Array2</code>.</p>\n<p><em><strong>See also all methods for <a href=\"ndarray/struct.ArrayBase.html\" title=\"struct ndarray::ArrayBase\"><code>ArrayBase</code></a></strong></em></p>\n</div><div class=\"impl-items\"><details class=\"toggle method-toggle\" open><summary><section id=\"method.push_row\" class=\"method\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#114-119\">source</a><h4 class=\"code-header\">pub fn <a href=\"ndarray/type.Array.html#tymethod.push_row\" class=\"fn\">push_row</a>(&amp;mut self, row: <a class=\"type\" href=\"ndarray/type.ArrayView.html\" title=\"type ndarray::ArrayView\">ArrayView</a>&lt;'_, A, <a class=\"type\" href=\"ndarray/type.Ix1.html\" title=\"type ndarray::Ix1\">Ix1</a>&gt;) -&gt; <a class=\"enum\" href=\"https://doc.rust-lang.org/1.77.0/core/result/enum.Result.html\" title=\"enum core::result::Result\">Result</a>&lt;<a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.unit.html\">()</a>, <a class=\"struct\" href=\"ndarray/struct.ShapeError.html\" title=\"struct ndarray::ShapeError\">ShapeError</a>&gt;<div class=\"where\">where\n    A: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/clone/trait.Clone.html\" title=\"trait core::clone::Clone\">Clone</a>,</div></h4></section></summary><div class=\"docblock\"><p>Append a row to an array</p>\n<p>The elements from <code>row</code> are cloned and added as a new row in the array.</p>\n<p><em><strong>Errors</strong></em> with a shape error if the length of the row does not match the length of the\nrows in the array.</p>\n<p>The memory layout of the <code>self</code> array matters for ensuring that the append is efficient.\nAppending automatically changes memory layout of the array so that it is appended to\nalong the “growing axis”. However, if the memory layout needs adjusting, the array must\nreallocate and move memory.</p>\n<p>The operation leaves the existing data in place and is most efficent if one of these is\ntrue:</p>\n<ul>\n<li>The axis being appended to is the longest stride axis, i.e the array is in row major\n(“C”) layout.</li>\n<li>The array has 0 or 1 rows (It is converted to row major)</li>\n</ul>\n<p>Ensure appending is efficient by, for example, appending to an empty array and then always\npushing/appending along the same axis. For pushing rows, ndarray’s default layout (C order)\nis efficient.</p>\n<p>When repeatedly appending to a single axis, the amortized average complexity of each\nappend is O(m), where <em>m</em> is the length of the row.</p>\n\n<div class=\"example-wrap\"><pre class=\"rust rust-example-rendered\"><code><span class=\"kw\">use </span>ndarray::{Array, ArrayView, array};\n\n<span class=\"comment\">// create an empty array and append\n</span><span class=\"kw\">let </span><span class=\"kw-2\">mut </span>a = Array::zeros((<span class=\"number\">0</span>, <span class=\"number\">4</span>));\na.push_row(ArrayView::from(<span class=\"kw-2\">&amp;</span>[ <span class=\"number\">1.</span>,  <span class=\"number\">2.</span>,  <span class=\"number\">3.</span>,  <span class=\"number\">4.</span>])).unwrap();\na.push_row(ArrayView::from(<span class=\"kw-2\">&amp;</span>[-<span class=\"number\">1.</span>, -<span class=\"number\">2.</span>, -<span class=\"number\">3.</span>, -<span class=\"number\">4.</span>])).unwrap();\n\n<span class=\"macro\">assert_eq!</span>(\n    a,\n    <span class=\"macro\">array!</span>[[ <span class=\"number\">1.</span>,  <span class=\"number\">2.</span>,  <span class=\"number\">3.</span>,  <span class=\"number\">4.</span>],\n           [-<span class=\"number\">1.</span>, -<span class=\"number\">2.</span>, -<span class=\"number\">3.</span>, -<span class=\"number\">4.</span>]]);</code></pre></div>\n</div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.push_column\" class=\"method\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#160-165\">source</a><h4 class=\"code-header\">pub fn <a href=\"ndarray/type.Array.html#tymethod.push_column\" class=\"fn\">push_column</a>(\n    &amp;mut self,\n    column: <a class=\"type\" href=\"ndarray/type.ArrayView.html\" title=\"type ndarray::ArrayView\">ArrayView</a>&lt;'_, A, <a class=\"type\" href=\"ndarray/type.Ix1.html\" title=\"type ndarray::Ix1\">Ix1</a>&gt;\n) -&gt; <a class=\"enum\" href=\"https://doc.rust-lang.org/1.77.0/core/result/enum.Result.html\" title=\"enum core::result::Result\">Result</a>&lt;<a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.unit.html\">()</a>, <a class=\"struct\" href=\"ndarray/struct.ShapeError.html\" title=\"struct ndarray::ShapeError\">ShapeError</a>&gt;<div class=\"where\">where\n    A: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/clone/trait.Clone.html\" title=\"trait core::clone::Clone\">Clone</a>,</div></h4></section></summary><div class=\"docblock\"><p>Append a column to an array</p>\n<p>The elements from <code>column</code> are cloned and added as a new column in the array.</p>\n<p><em><strong>Errors</strong></em> with a shape error if the length of the column does not match the length of\nthe columns in the array.</p>\n<p>The memory layout of the <code>self</code> array matters for ensuring that the append is efficient.\nAppending automatically changes memory layout of the array so that it is appended to\nalong the “growing axis”. However, if the memory layout needs adjusting, the array must\nreallocate and move memory.</p>\n<p>The operation leaves the existing data in place and is most efficent if one of these is\ntrue:</p>\n<ul>\n<li>The axis being appended to is the longest stride axis, i.e the array is in column major\n(“F”) layout.</li>\n<li>The array has 0 or 1 columns (It is converted to column major)</li>\n</ul>\n<p>Ensure appending is efficient by, for example, appending to an empty array and then always\npushing/appending along the same axis. For pushing columns, column major layout (F order)\nis efficient.</p>\n<p>When repeatedly appending to a single axis, the amortized average complexity of each append\nis O(m), where <em>m</em> is the length of the column.</p>\n\n<div class=\"example-wrap\"><pre class=\"rust rust-example-rendered\"><code><span class=\"kw\">use </span>ndarray::{Array, ArrayView, array};\n\n<span class=\"comment\">// create an empty array and append\n</span><span class=\"kw\">let </span><span class=\"kw-2\">mut </span>a = Array::zeros((<span class=\"number\">2</span>, <span class=\"number\">0</span>));\na.push_column(ArrayView::from(<span class=\"kw-2\">&amp;</span>[<span class=\"number\">1.</span>, <span class=\"number\">2.</span>])).unwrap();\na.push_column(ArrayView::from(<span class=\"kw-2\">&amp;</span>[-<span class=\"number\">1.</span>, -<span class=\"number\">2.</span>])).unwrap();\n\n<span class=\"macro\">assert_eq!</span>(\n    a,\n    <span class=\"macro\">array!</span>[[<span class=\"number\">1.</span>, -<span class=\"number\">1.</span>],\n           [<span class=\"number\">2.</span>, -<span class=\"number\">2.</span>]]);</code></pre></div>\n</div></details></div></details>",0,"ndarray::aliases::Array2"],["<details class=\"toggle implementors-toggle\" open><summary><section id=\"impl-ArrayBase%3COwnedRepr%3CA%3E,+D%3E\" class=\"impl\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#168-660\">source</a><a href=\"#impl-ArrayBase%3COwnedRepr%3CA%3E,+D%3E\" class=\"anchor\">§</a><h3 class=\"code-header\">impl&lt;A, D&gt; <a class=\"type\" href=\"ndarray/type.Array.html\" title=\"type ndarray::Array\">Array</a>&lt;A, D&gt;<div class=\"where\">where\n    D: <a class=\"trait\" href=\"ndarray/trait.Dimension.html\" title=\"trait ndarray::Dimension\">Dimension</a>,</div></h3></section></summary><div class=\"impl-items\"><details class=\"toggle method-toggle\" open><summary><section id=\"method.move_into\" class=\"method\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#189-205\">source</a><h4 class=\"code-header\">pub fn <a href=\"ndarray/type.Array.html#tymethod.move_into\" class=\"fn\">move_into</a>&lt;'a, AM&gt;(self, new_array: AM)<div class=\"where\">where\n    AM: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/convert/trait.Into.html\" title=\"trait core::convert::Into\">Into</a>&lt;<a class=\"type\" href=\"ndarray/type.ArrayViewMut.html\" title=\"type ndarray::ArrayViewMut\">ArrayViewMut</a>&lt;'a, A, D&gt;&gt;,\n    A: 'a,</div></h4></section></summary><div class=\"docblock\"><p>Move all elements from self into <code>new_array</code>, which must be of the same shape but\ncan have a different memory layout. The destination is overwritten completely.</p>\n<p>The destination should be a mut reference to an array or an <code>ArrayViewMut</code> with\n<code>A</code> elements.</p>\n<p><em><strong>Panics</strong></em> if the shapes don’t agree.</p>\n<h6 id=\"example\"><a class=\"doc-anchor\" href=\"#example\">§</a>Example</h6>\n<div class=\"example-wrap\"><pre class=\"rust rust-example-rendered\"><code><span class=\"kw\">use </span>ndarray::Array;\n\n<span class=\"comment\">// Usage example of move_into in safe code\n</span><span class=\"kw\">let </span><span class=\"kw-2\">mut </span>a = Array::default((<span class=\"number\">10</span>, <span class=\"number\">10</span>));\n<span class=\"kw\">let </span>b = Array::from_shape_fn((<span class=\"number\">10</span>, <span class=\"number\">10</span>), |(i, j)| (i + j).to_string());\nb.move_into(<span class=\"kw-2\">&amp;mut </span>a);</code></pre></div>\n</div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.move_into_uninit\" class=\"method\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#242-249\">source</a><h4 class=\"code-header\">pub fn <a href=\"ndarray/type.Array.html#tymethod.move_into_uninit\" class=\"fn\">move_into_uninit</a>&lt;'a, AM&gt;(self, new_array: AM)<div class=\"where\">where\n    AM: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/convert/trait.Into.html\" title=\"trait core::convert::Into\">Into</a>&lt;<a class=\"type\" href=\"ndarray/type.ArrayViewMut.html\" title=\"type ndarray::ArrayViewMut\">ArrayViewMut</a>&lt;'a, <a class=\"union\" href=\"https://doc.rust-lang.org/1.77.0/core/mem/maybe_uninit/union.MaybeUninit.html\" title=\"union core::mem::maybe_uninit::MaybeUninit\">MaybeUninit</a>&lt;A&gt;, D&gt;&gt;,\n    A: 'a,</div></h4></section></summary><div class=\"docblock\"><p>Move all elements from self into <code>new_array</code>, which must be of the same shape but\ncan have a different memory layout. The destination is overwritten completely.</p>\n<p>The destination should be a mut reference to an array or an <code>ArrayViewMut</code> with\n<code>MaybeUninit&lt;A&gt;</code> elements (which are overwritten without dropping any existing value).</p>\n<p>Minor implementation note: Owned arrays like <code>self</code> may be sliced in place and own elements\nthat are not part of their active view; these are dropped at the end of this function,\nafter all elements in the “active view” are moved into <code>new_array</code>. If there is a panic in\ndrop of any such element, other elements may be leaked.</p>\n<p><em><strong>Panics</strong></em> if the shapes don’t agree.</p>\n<h6 id=\"example-1\"><a class=\"doc-anchor\" href=\"#example-1\">§</a>Example</h6>\n<div class=\"example-wrap\"><pre class=\"rust rust-example-rendered\"><code><span class=\"kw\">use </span>ndarray::Array;\n\n<span class=\"kw\">let </span>a = Array::from_iter(<span class=\"number\">0</span>..<span class=\"number\">100</span>).into_shape((<span class=\"number\">10</span>, <span class=\"number\">10</span>)).unwrap();\n<span class=\"kw\">let </span><span class=\"kw-2\">mut </span>b = Array::uninit((<span class=\"number\">10</span>, <span class=\"number\">10</span>));\na.move_into_uninit(<span class=\"kw-2\">&amp;mut </span>b);\n<span class=\"kw\">unsafe </span>{\n    <span class=\"comment\">// we can now promise we have fully initialized `b`.\n    </span><span class=\"kw\">let </span>b = b.assume_init();\n}</code></pre></div>\n</div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.push\" class=\"method\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#404-412\">source</a><h4 class=\"code-header\">pub fn <a href=\"ndarray/type.Array.html#tymethod.push\" class=\"fn\">push</a>(\n    &amp;mut self,\n    axis: <a class=\"struct\" href=\"ndarray/struct.Axis.html\" title=\"struct ndarray::Axis\">Axis</a>,\n    array: <a class=\"type\" href=\"ndarray/type.ArrayView.html\" title=\"type ndarray::ArrayView\">ArrayView</a>&lt;'_, A, D::<a class=\"associatedtype\" href=\"ndarray/trait.Dimension.html#associatedtype.Smaller\" title=\"type ndarray::Dimension::Smaller\">Smaller</a>&gt;\n) -&gt; <a class=\"enum\" href=\"https://doc.rust-lang.org/1.77.0/core/result/enum.Result.html\" title=\"enum core::result::Result\">Result</a>&lt;<a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.unit.html\">()</a>, <a class=\"struct\" href=\"ndarray/struct.ShapeError.html\" title=\"struct ndarray::ShapeError\">ShapeError</a>&gt;<div class=\"where\">where\n    A: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/clone/trait.Clone.html\" title=\"trait core::clone::Clone\">Clone</a>,\n    D: <a class=\"trait\" href=\"ndarray/trait.RemoveAxis.html\" title=\"trait ndarray::RemoveAxis\">RemoveAxis</a>,</div></h4></section></summary><div class=\"docblock\"><p>Append an array to the array along an axis.</p>\n<p>The elements of <code>array</code> are cloned and extend the axis <code>axis</code> in the present array;\n<code>self</code> will grow in size by 1 along <code>axis</code>.</p>\n<p>Append to the array, where the array being pushed to the array has one dimension less than\nthe <code>self</code> array. This method is equivalent to <a href=\"ndarray/struct.ArrayBase.html#method.append\" title=\"method ndarray::ArrayBase::append\">append</a> in this way:\n<code>self.append(axis, array.insert_axis(axis))</code>.</p>\n<p><em><strong>Errors</strong></em> with a shape error if the shape of self does not match the array-to-append;\nall axes <em>except</em> the axis along which it being appended matter for this check:\nthe shape of <code>self</code> with <code>axis</code> removed must be the same as the shape of <code>array</code>.</p>\n<p>The memory layout of the <code>self</code> array matters for ensuring that the append is efficient.\nAppending automatically changes memory layout of the array so that it is appended to\nalong the “growing axis”. However, if the memory layout needs adjusting, the array must\nreallocate and move memory.</p>\n<p>The operation leaves the existing data in place and is most efficent if <code>axis</code> is a\n“growing axis” for the array, i.e. one of these is true:</p>\n<ul>\n<li>The axis is the longest stride axis, for example the 0th axis in a C-layout or the\n<em>n-1</em>th axis in an F-layout array.</li>\n<li>The axis has length 0 or 1 (It is converted to the new growing axis)</li>\n</ul>\n<p>Ensure appending is efficient by for example starting from an empty array and/or always\nappending to an array along the same axis.</p>\n<p>The amortized average complexity of the append, when appending along its growing axis, is\nO(<em>m</em>) where <em>m</em> is the number of individual elements to append.</p>\n<p>The memory layout of the argument <code>array</code> does not matter to the same extent.</p>\n\n<div class=\"example-wrap\"><pre class=\"rust rust-example-rendered\"><code><span class=\"kw\">use </span>ndarray::{Array, ArrayView, array, Axis};\n\n<span class=\"comment\">// create an empty array and push rows to it\n</span><span class=\"kw\">let </span><span class=\"kw-2\">mut </span>a = Array::zeros((<span class=\"number\">0</span>, <span class=\"number\">4</span>));\n<span class=\"kw\">let </span>ones  = ArrayView::from(<span class=\"kw-2\">&amp;</span>[<span class=\"number\">1.</span>; <span class=\"number\">4</span>]);\n<span class=\"kw\">let </span>zeros = ArrayView::from(<span class=\"kw-2\">&amp;</span>[<span class=\"number\">0.</span>; <span class=\"number\">4</span>]);\na.push(Axis(<span class=\"number\">0</span>), ones).unwrap();\na.push(Axis(<span class=\"number\">0</span>), zeros).unwrap();\na.push(Axis(<span class=\"number\">0</span>), ones).unwrap();\n\n<span class=\"macro\">assert_eq!</span>(\n    a,\n    <span class=\"macro\">array!</span>[[<span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>],\n           [<span class=\"number\">0.</span>, <span class=\"number\">0.</span>, <span class=\"number\">0.</span>, <span class=\"number\">0.</span>],\n           [<span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>]]);</code></pre></div>\n</div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.append\" class=\"method\"><a class=\"src rightside\" href=\"src/ndarray/impl_owned_array.rs.html#465-659\">source</a><h4 class=\"code-header\">pub fn <a href=\"ndarray/type.Array.html#tymethod.append\" class=\"fn\">append</a>(\n    &amp;mut self,\n    axis: <a class=\"struct\" href=\"ndarray/struct.Axis.html\" title=\"struct ndarray::Axis\">Axis</a>,\n    array: <a class=\"type\" href=\"ndarray/type.ArrayView.html\" title=\"type ndarray::ArrayView\">ArrayView</a>&lt;'_, A, D&gt;\n) -&gt; <a class=\"enum\" href=\"https://doc.rust-lang.org/1.77.0/core/result/enum.Result.html\" title=\"enum core::result::Result\">Result</a>&lt;<a class=\"primitive\" href=\"https://doc.rust-lang.org/1.77.0/std/primitive.unit.html\">()</a>, <a class=\"struct\" href=\"ndarray/struct.ShapeError.html\" title=\"struct ndarray::ShapeError\">ShapeError</a>&gt;<div class=\"where\">where\n    A: <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/clone/trait.Clone.html\" title=\"trait core::clone::Clone\">Clone</a>,\n    D: <a class=\"trait\" href=\"ndarray/trait.RemoveAxis.html\" title=\"trait ndarray::RemoveAxis\">RemoveAxis</a>,</div></h4></section></summary><div class=\"docblock\"><p>Append an array to the array along an axis.</p>\n<p>The elements of <code>array</code> are cloned and extend the axis <code>axis</code> in the present array;\n<code>self</code> will grow in size by <code>array.len_of(axis)</code> along <code>axis</code>.</p>\n<p><em><strong>Errors</strong></em> with a shape error if the shape of self does not match the array-to-append;\nall axes <em>except</em> the axis along which it being appended matter for this check:\nthe shape of <code>self</code> with <code>axis</code> removed must be the same as the shape of <code>array</code> with\n<code>axis</code> removed.</p>\n<p>The memory layout of the <code>self</code> array matters for ensuring that the append is efficient.\nAppending automatically changes memory layout of the array so that it is appended to\nalong the “growing axis”. However, if the memory layout needs adjusting, the array must\nreallocate and move memory.</p>\n<p>The operation leaves the existing data in place and is most efficent if <code>axis</code> is a\n“growing axis” for the array, i.e. one of these is true:</p>\n<ul>\n<li>The axis is the longest stride axis, for example the 0th axis in a C-layout or the\n<em>n-1</em>th axis in an F-layout array.</li>\n<li>The axis has length 0 or 1 (It is converted to the new growing axis)</li>\n</ul>\n<p>Ensure appending is efficient by for example starting from an empty array and/or always\nappending to an array along the same axis.</p>\n<p>The amortized average complexity of the append, when appending along its growing axis, is\nO(<em>m</em>) where <em>m</em> is the number of individual elements to append.</p>\n<p>The memory layout of the argument <code>array</code> does not matter to the same extent.</p>\n\n<div class=\"example-wrap\"><pre class=\"rust rust-example-rendered\"><code><span class=\"kw\">use </span>ndarray::{Array, ArrayView, array, Axis};\n\n<span class=\"comment\">// create an empty array and append two rows at a time\n</span><span class=\"kw\">let </span><span class=\"kw-2\">mut </span>a = Array::zeros((<span class=\"number\">0</span>, <span class=\"number\">4</span>));\n<span class=\"kw\">let </span>ones  = ArrayView::from(<span class=\"kw-2\">&amp;</span>[<span class=\"number\">1.</span>; <span class=\"number\">8</span>]).into_shape((<span class=\"number\">2</span>, <span class=\"number\">4</span>)).unwrap();\n<span class=\"kw\">let </span>zeros = ArrayView::from(<span class=\"kw-2\">&amp;</span>[<span class=\"number\">0.</span>; <span class=\"number\">8</span>]).into_shape((<span class=\"number\">2</span>, <span class=\"number\">4</span>)).unwrap();\na.append(Axis(<span class=\"number\">0</span>), ones).unwrap();\na.append(Axis(<span class=\"number\">0</span>), zeros).unwrap();\na.append(Axis(<span class=\"number\">0</span>), ones).unwrap();\n\n<span class=\"macro\">assert_eq!</span>(\n    a,\n    <span class=\"macro\">array!</span>[[<span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>],\n           [<span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>],\n           [<span class=\"number\">0.</span>, <span class=\"number\">0.</span>, <span class=\"number\">0.</span>, <span class=\"number\">0.</span>],\n           [<span class=\"number\">0.</span>, <span class=\"number\">0.</span>, <span class=\"number\">0.</span>, <span class=\"number\">0.</span>],\n           [<span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>],\n           [<span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>, <span class=\"number\">1.</span>]]);</code></pre></div>\n</div></details></div></details>",0,"ndarray::aliases::Array0","ndarray::aliases::Array1","ndarray::aliases::Array2","ndarray::aliases::Array3","ndarray::aliases::Array4","ndarray::aliases::Array5","ndarray::aliases::Array6","ndarray::aliases::ArrayD"],["<details class=\"toggle implementors-toggle\" open><summary><section id=\"impl-IntoIterator-for-ArrayBase%3COwnedRepr%3CA%3E,+D%3E\" class=\"impl\"><a class=\"src rightside\" href=\"src/ndarray/iterators/into_iter.rs.html#100-110\">source</a><a href=\"#impl-IntoIterator-for-ArrayBase%3COwnedRepr%3CA%3E,+D%3E\" class=\"anchor\">§</a><h3 class=\"code-header\">impl&lt;A, D&gt; <a class=\"trait\" href=\"https://doc.rust-lang.org/1.77.0/core/iter/traits/collect/trait.IntoIterator.html\" title=\"trait core::iter::traits::collect::IntoIterator\">IntoIterator</a> for <a class=\"type\" href=\"ndarray/type.Array.html\" title=\"type ndarray::Array\">Array</a>&lt;A, D&gt;<div class=\"where\">where\n    D: <a class=\"trait\" href=\"ndarray/trait.Dimension.html\" title=\"trait ndarray::Dimension\">Dimension</a>,</div></h3></section></summary><div class=\"impl-items\"><details class=\"toggle\" open><summary><section id=\"associatedtype.Item\" class=\"associatedtype trait-impl\"><a href=\"#associatedtype.Item\" class=\"anchor\">§</a><h4 class=\"code-header\">type <a href=\"https://doc.rust-lang.org/1.77.0/core/iter/traits/collect/trait.IntoIterator.html#associatedtype.Item\" class=\"associatedtype\">Item</a> = A</h4></section></summary><div class='docblock'>The type of the elements being iterated over.</div></details><details class=\"toggle\" open><summary><section id=\"associatedtype.IntoIter\" class=\"associatedtype trait-impl\"><a href=\"#associatedtype.IntoIter\" class=\"anchor\">§</a><h4 class=\"code-header\">type <a href=\"https://doc.rust-lang.org/1.77.0/core/iter/traits/collect/trait.IntoIterator.html#associatedtype.IntoIter\" class=\"associatedtype\">IntoIter</a> = IntoIter&lt;A, D&gt;</h4></section></summary><div class='docblock'>Which kind of iterator are we turning this into?</div></details><details class=\"toggle method-toggle\" open><summary><section id=\"method.into_iter\" class=\"method trait-impl\"><a class=\"src rightside\" href=\"src/ndarray/iterators/into_iter.rs.html#107-109\">source</a><a href=\"#method.into_iter\" class=\"anchor\">§</a><h4 class=\"code-header\">fn <a href=\"https://doc.rust-lang.org/1.77.0/core/iter/traits/collect/trait.IntoIterator.html#tymethod.into_iter\" class=\"fn\">into_iter</a>(self) -&gt; Self::<a class=\"associatedtype\" href=\"https://doc.rust-lang.org/1.77.0/core/iter/traits/collect/trait.IntoIterator.html#associatedtype.IntoIter\" title=\"type core::iter::traits::collect::IntoIterator::IntoIter\">IntoIter</a></h4></section></summary><div class='docblock'>Creates an iterator from a value. <a href=\"https://doc.rust-lang.org/1.77.0/core/iter/traits/collect/trait.IntoIterator.html#tymethod.into_iter\">Read more</a></div></details></div></details>","IntoIterator","ndarray::aliases::Array0","ndarray::aliases::Array1","ndarray::aliases::Array2","ndarray::aliases::Array3","ndarray::aliases::Array4","ndarray::aliases::Array5","ndarray::aliases::Array6","ndarray::aliases::ArrayD"]]
};if (window.register_type_impls) {window.register_type_impls(type_impls);} else {window.pending_type_impls = type_impls;}})()