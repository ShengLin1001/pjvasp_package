"""Publisher-site knowledge for locating an article's official PDF.

Everything here is browser-agnostic: URL rules, page-state classification,
candidate filtering, and the JavaScript snippets a driver injects. The
browser driving itself stays in the download orchestrator, so the same
knowledge can back a Playwright, CDP, or CloakBrowser front end.

Functions:
    get_page_state: Classify a publisher page from its HTML and final URL.
    get_publisher_pdf_url: Derive a PDF URL from an article URL by host rule.
    check_supplement_url: Detect supplementary-material links.
    check_article_url: Detect links belonging to the requested article.
    filter_pdf_candidates: Drop supplements, cited papers, and foreign hosts.

Constants:
    JS_IS_PDF_DOCUMENT: Detect a tab currently displaying a PDF.
"""

from __future__ import annotations

import re
from urllib.parse import urlsplit


### page state, to here ###

def get_page_state(html: str, url: str) -> str:
    """Classify what a publisher page is currently showing.

    Args:
        html: Full page HTML; only the first 200k characters are scanned.
        url: The tab's final URL after redirects.

    Returns:
        One of ``pdf_ready``, ``cloudflare``, ``captcha``, ``article``,
        ``paywall``, ``institution_login``, or ``unknown``.
    """
    low = (html or "")[:200000].lower()
    if re.search(r"\.pdf(\?|#|$)", url or "", re.I):
        return "pdf_ready"
    # Cloudflare ships several interstitial wordings; match the markup too, so
    # a reworded challenge page is not mistaken for a loaded article.
    if ("just a moment" in low or "cf-turnstile" in low or "cf_chl_opt" in low
            or "challenges.cloudflare.com" in low or "verifying you are human" in low
            or "request verification" in low
            or "review the security of your connection" in low
            or "enable javascript and cookies to continue" in low):
        return "cloudflare"
    if "hcaptcha" in low or "g-recaptcha" in low or "recaptcha/api" in low:
        return "captcha"
    if "citation_pdf_url" in low or re.search(r">\s*(download|view)\s*(full[\s-]*text\s*)?pdf", low):
        return "article"
    if "get access" in low or "purchase pdf" in low or "buy article" in low:
        return "paywall"
    if "sign in via your institution" in low or "access through your institution" in low:
        return "institution_login"
    return "unknown"


### PDF URL rules, to here ###

# Nearly every publisher emits <meta name="citation_pdf_url"> for Google
# Scholar, so that tag replaces a table of per-host rules. Only these hosts
# need a hint: they either omit the tag, or advertise a URL that does not
# actually serve the file.
LHOST_PDF_RULES = (
    # ScienceDirect hides the PDF behind a viewer page; the asset is pii + /pdfft.
    (r"^(https://[^/]*sciencedirect\.com/science/article/pii/[A-Z0-9]+).*",
     r"\1/pdfft?isDTMRedir=true&download=true"),
    # APS advertises link.aps.org/pdf/..., which only bounces back to the
    # abstract. The journal host serves the real file.
    (r"^(https://journals\.aps\.org/[^/]+)/abstract/(10\..+)$", r"\1/pdf/\2"),
    (r"^(https://dl\.acm\.org)/doi/(?:abs/|full/)?(10\..+)$", r"\1/doi/pdf/\2"),
    # Wiley's /doi/pdf/ is a viewer wrapper; pdfdirect is the file itself.
    (r"^(https://[^/]*onlinelibrary\.wiley\.com)/doi/(?:abs/|full/|epdf/|pdf/)?(10\..+?)(?:\?.*)?$",
     r"\1/doi/pdfdirect/\2?download=true"),
    (r"^(https://pubs\.acs\.org)/doi/(?:abs/|full/|epdf/)?(10\..+?)(?:\?.*)?$", r"\1/doi/pdf/\2"),
    # Annual Reviews publishes no PDF link in the article DOM at all; its
    # Atypon crawler endpoint serves the file.
    (r"^(https://[^/]*annualreviews\.org/content/journals/10\.[^?]+?)(?:\?.*)?$",
     r"\1?crawler=true&mimetype=application/pdf"),
)


def get_publisher_pdf_url(url: str) -> str | None:
    """Derive a PDF URL from an article URL, for hosts whose metadata lies.

    Args:
        url: The article page URL after redirects.

    Returns:
        A PDF URL, or ``None`` when no host rule applies.
    """
    for pattern, replacement in LHOST_PDF_RULES:
        derived, count = re.subn(pattern, replacement, url or "")
        if count:
            return derived
    return None


### candidate filtering, to here ###

# Supplementary material, not the article. Publishers link it right beside the
# real PDF - ACS as "article-supplement", Springer under /esm/, OUP as
# "_supplemental_file" - and a plain href scan grabs it first.
SUPPLEMENT_PATTERN = re.compile(
    r"article-supplement|/esm/|_si_\d|suppl(_file|emental|ementary)|supporting-information",
    re.I)


def check_supplement_url(url: str) -> bool:
    """Return whether a URL points at supplementary material."""
    return bool(SUPPLEMENT_PATTERN.search(url or ""))


def check_article_url(url: str, doi: str, page_url: str) -> bool:
    """Return whether a URL belongs to the requested article.

    An article page links the PDFs of everything it cites, so a bare href
    scan returns neighbouring papers - and, on AIP, a product manual. Match
    on the DOI suffix or on the publisher's own article id in the page URL
    (Elsevier's PII).

    Args:
        url: Candidate PDF URL.
        doi: Requested DOI, e.g. ``10.1038/nature12373``.
        page_url: URL of the article page the candidate was found on.
    """
    low = (url or "").lower()
    if doi and doi.split("/", 1)[-1].lower() in low:
        return True
    match = re.search(r"/pii/([A-Z0-9]+)", page_url or "", re.I)
    return bool(match and match.group(1).lower() in low)


def filter_pdf_candidates(lurl, doi: str, page_url: str) -> list[str]:
    """Drop supplements, cited papers, and foreign hosts from PDF candidates.

    Args:
        lurl: Candidate PDF URLs, best first.
        doi: Requested DOI.
        page_url: URL of the article page.

    Returns:
        Candidates worth trying, order preserved. Falls back to same-host
        candidates when none carry the article's own id, because some
        publishers (MDPI, ACS) use ids unrelated to the DOI - but never to a
        foreign host, which is only ever a cited paper or a marketing PDF.
    """
    lkept = []
    for url in lurl:
        if url and url not in lkept and not check_supplement_url(url):
            lkept.append(url)
    lown = [url for url in lkept if check_article_url(url, doi, page_url)]
    if lown:
        return lown
    host = urlsplit(page_url or "").netloc
    return [url for url in lkept if urlsplit(url).netloc == host]


### injected JavaScript, to here ###

# Collect PDF candidates from the article DOM, best first.
JS_FIND_PDF_URLS = r"""
(() => {
  const abs = (h) => { try { return new URL(h, location.href).href; } catch (e) { return null; } };
  const out = [];
  const push = (u) => { if (u && !out.includes(u)) out.push(u); };
  const meta = document.querySelector('meta[name="citation_pdf_url"]');
  if (meta) push(abs(meta.content));
  // IEEE links stamp.jsp, but that is only a wrapper whose iframe holds the
  // real file at stampPDF/getPDF.jsp. Offer the inner URL first.
  const stamp = document.querySelector('a[href*="stamp.jsp"],iframe[src*="stamp.jsp"]');
  if (stamp) {
    const s = abs(stamp.getAttribute('href') || stamp.getAttribute('src'));
    if (s) { push(s.replace('/stamp/stamp.jsp', '/stampPDF/getPDF.jsp')); push(s); }
  }
  const links = [...document.querySelectorAll('a[href]')];
  const txt = (e) => ((e.innerText||'') + ' ' + (e.getAttribute('aria-label')||''))
      .replace(/\s+/g, ' ').trim().toLowerCase();
  // A labelled button first: "Standard PDF" beats the enhanced-reader link.
  for (const a of links)
    if (/^(standard pdf|download pdf|view pdf|full[- ]?text pdf|article pdf|pdf \(\d)/.test(txt(a)))
      push(abs(a.getAttribute('href')));
  for (const a of links) {
    const h = a.getAttribute('href');
    if (h && (/\.pdf(\?|#|$)/i.test(h) || /\/(pdf|pdfdirect|pdfft|printable)\b/i.test(h)))
      push(abs(h));
  }
  return out.slice(0, 5);
})()
"""

# Detect a tab that is *showing* a PDF, however it got there. A URL suffix is
# not enough: plenty of publishers serve the file from an extension-less path
# (PNAS /doi/pdf/, IEEE getPDF.jsp, Annual Reviews ?mimetype=application/pdf).
JS_IS_PDF_DOCUMENT = r"""
(() => document.contentType === 'application/pdf'
     || !!document.querySelector(
          'embed[type="application/pdf"],object[type="application/pdf"]'))()
"""

# Cookie-consent overlays sit on top of the article (AIP, IEEE) and block both
# clicks and the PDF request.
JS_ACCEPT_CONSENT = r"""
(() => {
  const re = /^(accept all|accept cookies|accept & close|i accept|accept|agree|allow all|got it|同意|接受)/i;
  const t = [...document.querySelectorAll('button,a,[role=button]')]
      .find(e => re.test((e.innerText||'').replace(/\s+/g,' ').trim()) && e.offsetParent !== null);
  if (t) { t.click(); return true; }
  return false;
})()
"""

# Open the publisher's institutional-access entry point. Elsevier and Cell
# Press hide it behind "Access through your institution"; Atypon hosts say
# "Institutional Login"; some show only a bare "Sign in".
JS_OPEN_INSTITUTION = r"""
(() => {
  const re = /access through your institution|institutional (log ?in|access|sign)|sign in via your institution|log in via (your )?institution|通过您的机构/i;
  const els = [...document.querySelectorAll('a,button,[role=button]')];
  const t = els.find(e => re.test((e.innerText||'') + ' ' + (e.getAttribute('aria-label')||'')));
  if (t) { t.click(); return 'inst'; }
  const g = els.find(e => /^(get access|access|sign in)$/i.test((e.innerText||'').trim()));
  if (g) { g.click(); return 'access'; }
  return null;
})()
"""

# The entry point opens a typeahead of institutions: type, then pick.
JS_TYPE_INSTITUTION = r"""
(name) => {
  const box = document.querySelector(
    'input[type=search],input[placeholder*="nstitution" i],input[placeholder*="rgani" i],input[id*="nstitution" i]');
  if (!box) return null;
  const set = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value').set;
  set.call(box, name);
  box.dispatchEvent(new Event('input', {bubbles: true}));
  box.dispatchEvent(new Event('change', {bubbles: true}));
  return 'typed';
}
"""

JS_PICK_INSTITUTION = r"""
(name) => {
  const els = [...document.querySelectorAll('a,button,li,[role=option],[role=button]')];
  const t = els.find(e => (e.innerText||'').trim().toLowerCase().includes(name.toLowerCase()));
  if (t) { t.click(); return 'picked'; }
  return null;
}
"""

# Read a fetched PDF back out of the page in base64 chunks. Each chunk must be
# decoded on its own: every chunk carries its own padding, so concatenating the
# *encoded* strings corrupts the file.
JS_FETCH_PDF = """
async (url) => {
  const r = await fetch(url, {credentials: 'include', redirect: 'follow'});
  if (!r.ok) throw new Error('http ' + r.status);
  window.__buf = new Uint8Array(await r.arrayBuffer());
  return window.__buf.length;
}
"""

JS_READ_CHUNK = """
([off, len]) => {
  const u = window.__buf.subarray(off, off + len);
  let s = '';
  for (let i = 0; i < u.length; i += 0x8000)
    s += String.fromCharCode.apply(null, u.subarray(i, i + 0x8000));
  return btoa(s);
}
"""
