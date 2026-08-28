"""Tests for reusable literature metadata and filename helpers."""

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


if __name__ == "__main__":
    unittest.main()
