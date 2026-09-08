"""Tests for reusable literature metadata, filename, and publisher helpers."""

import tempfile
import unittest
from pathlib import Path

from mymetal.academic.search.literature_download import (
    check_journal_metadata,
    generate_pdf_filename,
    get_journal_abbreviation,
    is_complete_pdf,
    normalize_doi,
    parse_dois,
)
from mymetal.academic.search.publisher_pdf import (
    check_article_url,
    check_supplement_url,
    filter_pdf_candidates,
    get_page_state,
    get_publisher_pdf_url,
)


class TestLiterature(unittest.TestCase):
    def test_known_journal_aliases(self):
        self.assertEqual(get_journal_abbreviation("Physical Review B"), "PRB")
        self.assertEqual(get_journal_abbreviation("PRB"), "PRB")
        self.assertEqual(get_journal_abbreviation("Applied Physics Letters"), "APL")
        self.assertEqual(get_journal_abbreviation("Acta Materialia"), "ACTA-MATER")

    def test_generates_requested_filename(self):
        dict_metadata = {
            "type": "journal-article",
            "title": ["Effect of strain on the stacking fault energy of copper: A first-principles study"],
            "container-title": ["Physical Review B"],
            "published-print": {"date-parts": [[2013, 8, 26]]},
        }

        self.assertEqual(
            generate_pdf_filename(dict_metadata),
            "2013-PRB-Effect-of-st.pdf",
        )

    def test_keeps_first_ten_chinese_characters(self):
        dict_metadata = {
            "type": "journal-article",
            "title": ["材料计算模拟方法与应用进展"],
            "container-title": ["Acta Materialia"],
            "published-online": {"date-parts": [[2026]]},
        }

        self.assertEqual(
            generate_pdf_filename(dict_metadata),
            "2026-ACTA-MATER-材料计算模拟方法与应.pdf",
        )

    def test_rejects_preprint_snapshot_and_unknown_journal(self):
        self.assertEqual(
            check_journal_metadata({"type": "posted-content"}),
            "not a journal article",
        )
        self.assertEqual(
            check_journal_metadata({
                "type": "journal-article",
                "title": ["SnapShot: Mechanical Forces in Development I"],
                "container-title": ["Cell"],
            }),
            "excluded article type",
        )
        self.assertEqual(
            check_journal_metadata({
                "type": "journal-article",
                "title": ["A paper"],
                "container-title": ["Unknown Journal"],
            }),
            "journal abbreviation missing",
        )

    def test_parses_utf8_doi_list(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_file = Path(tmp) / "dois.txt"
            path_file.write_text(
                "\ufeffhttps://doi.org/10.1016/a\n# comment\n10.1016/a\n10.1103/b note\n",
                encoding="utf-8",
            )
            self.assertEqual(parse_dois(path_file), ["10.1016/a", "10.1103/b"])
            self.assertEqual(normalize_doi("doi: 10.1103/b"), "10.1103/b")

    def test_checks_complete_pdf_tail(self):
        with tempfile.TemporaryDirectory() as tmp:
            path_pdf = Path(tmp) / "paper.pdf"
            path_pdf.write_bytes(b"%PDF-" + b"x" * 5001 + b"%%EOF")
            self.assertTrue(is_complete_pdf(path_pdf))
            path_pdf.write_bytes(b"%PDF-" + b"x" * 5001)
            self.assertFalse(is_complete_pdf(path_pdf))


class TestPublisherPdf(unittest.TestCase):
    def test_classifies_publisher_pages(self):
        self.assertEqual(get_page_state("<html>Just a moment...</html>", "https://x/"),
                         "cloudflare")
        self.assertEqual(get_page_state("<div>Verifying you are human</div>", "https://x/"),
                         "cloudflare")
        self.assertEqual(get_page_state("<h1>Request Verification: In Progress</h1>",
                                        "https://x/"), "cloudflare")
        self.assertEqual(get_page_state("<html>hcaptcha</html>", "https://x/"), "captcha")
        self.assertEqual(get_page_state("<html/>", "https://x/a/1.pdf"), "pdf_ready")
        self.assertEqual(get_page_state('<meta name="citation_pdf_url" content="x">',
                                        "https://x/a"), "article")
        self.assertEqual(get_page_state("<html>Get access to this article</html>",
                                        "https://x/"), "paywall")

    def test_derives_pdf_url_for_hosts_whose_metadata_lies(self):
        self.assertEqual(
            get_publisher_pdf_url(
                "https://www.sciencedirect.com/science/article/pii/S0092867416303178"),
            "https://www.sciencedirect.com/science/article/pii/S0092867416303178"
            "/pdfft?isDTMRedir=true&download=true")
        self.assertEqual(
            get_publisher_pdf_url(
                "https://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.064104"),
            "https://journals.aps.org/prb/pdf/10.1103/PhysRevB.88.064104")
        self.assertEqual(
            get_publisher_pdf_url("https://dl.acm.org/doi/10.1145/3806644"),
            "https://dl.acm.org/doi/pdf/10.1145/3806644")
        self.assertEqual(
            get_publisher_pdf_url(
                "https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/adma.73337"),
            "https://advanced.onlinelibrary.wiley.com/doi/pdfdirect/10.1002/adma.73337"
            "?download=true")
        self.assertEqual(
            get_publisher_pdf_url("https://pubs.acs.org/doi/10.1021/jacs.5c20581"),
            "https://pubs.acs.org/doi/pdf/10.1021/jacs.5c20581")
        self.assertEqual(
            get_publisher_pdf_url(
                "https://www.annualreviews.org/content/journals/10.1146/annurev.psych.52.1.1"),
            "https://www.annualreviews.org/content/journals/10.1146/annurev.psych.52.1.1"
            "?crawler=true&mimetype=application/pdf")
        # A host that advertises a usable citation_pdf_url needs no rule.
        self.assertIsNone(
            get_publisher_pdf_url("https://www.nature.com/articles/nature12373"))

    def test_rejects_supplementary_material(self):
        for url in (
                "https://pubs.acs.org/jacsat/article-supplement/5098050/pdf/ja5c20581_si_001/",
                "https://media.springernature.com/springer-static/esm/art%3A10.1186/MOESM1.pdf",
                "https://oup.silverchair-cdn.com/x/gkag457_supplemental_file.pdf"):
            self.assertTrue(check_supplement_url(url), url)
        self.assertFalse(check_supplement_url(
            "https://link.springer.com/content/pdf/10.1186/1471-2105-15-135.pdf"))

    def test_keeps_only_this_articles_links(self):
        page = "https://www.sciencedirect.com/science/article/pii/S0927025618307924"
        doi = "10.1016/j.commatsci.2018.12.013"
        self.assertTrue(check_article_url(page + "/pdfft?isDTMRedir=true", doi, page))
        # A cited paper on the same host must not be mistaken for the article.
        self.assertFalse(check_article_url(
            "https://www.sciencedirect.com/science/article/pii/007964259290003P/pdf?md5=x",
            doi, page))
        self.assertEqual(
            filter_pdf_candidates(
                [page + "/pdfft?isDTMRedir=true",
                 "https://www.sciencedirect.com/science/article/pii/007964259290003P/pdf"],
                doi, page),
            [page + "/pdfft?isDTMRedir=true"])

    def test_falls_back_to_same_host_never_foreign(self):
        # MDPI ids are unrelated to the DOI, so the same-host fallback keeps
        # its link - but a marketing PDF on another host is always dropped.
        page = "https://www.mdpi.com/1424-8220/21/5/1705"
        self.assertEqual(
            filter_pdf_candidates(
                ["https://www.mdpi.com/1424-8220/21/5/1705/pdf?version=1",
                 "https://publishing.example.org/wp-content/uploads/Guide.pdf"],
                "10.3390/s21051705", page),
            ["https://www.mdpi.com/1424-8220/21/5/1705/pdf?version=1"])


if __name__ == "__main__":
    unittest.main()
