/* Draw-page math rendering: gui.js (vendored byte-identical from docs/gui/)
   typesets its ".arithmatex" blocks by calling window.MathJax.typesetPromise()
   -- the docs site ships MathJax, this webapp ships vendored KaTeX (strict CSP,
   no CDN). Bridge the one call gui.js makes onto KaTeX's auto-render. */
"use strict";
window.MathJax = {
  typesetPromise: function () {
    var root = document.getElementById("qlgui");
    if (root && window.renderMathInElement) {
      window.renderMathInElement(root, {
        delimiters: [
          { left: "\\[", right: "\\]", display: true },
          { left: "\\(", right: "\\)", display: false },
        ],
        throwOnError: false,
      });
    }
    return Promise.resolve();
  },
};
