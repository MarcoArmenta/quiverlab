# Vendored static assets

The strict CSP (`script-src 'self'`) forbids any runtime CDN, so KaTeX is
**self-hosted**: it is committed here and served by the app under
`/static/katex/`. No page ever hotlinks a CDN.

## KaTeX

- **Version:** 0.16.22
- **Upstream:** https://github.com/KaTeX/KaTeX/releases/tag/v0.16.22
  (npm package `katex@0.16.22`, `dist/`)
- **License:** MIT (KaTeX contributors)
- **Files committed** under `katex/`:
  - `katex.min.css` — stylesheet (references `fonts/` relative, same-origin)
  - `katex.min.js` — renderer
  - `auto-render.min.js` — the `contrib/auto-render` helper (`renderMathInElement`),
    invoked from `app.js` (CSP-safe, external)
  - `fonts/` — the 60 KaTeX web-font files (`.woff2`/`.woff`/`.ttf`)

### Re-vendoring at build/deploy time

To fetch a clean copy (e.g. in the Dockerfile image build) instead of relying on
the committed assets:

    mkdir -p webapp/static/katex
    curl -L https://github.com/KaTeX/KaTeX/releases/download/v0.16.22/katex.tar.gz \
      | tar xz --strip-components=1 -C webapp/static/katex \
      katex/katex.min.css katex/katex.min.js katex/contrib/auto-render.min.js katex/fonts

If KaTeX cannot be obtained, the pages still render: `app.js` guards the
`renderMathInElement` call, so math falls back to plain text rather than failing.

## app.js

The page client (form submit → `/api/compute`, big-job email flow, job polling,
KaTeX auto-render). External static file — never inlined — because the CSP
forbids inline `<script>`.
