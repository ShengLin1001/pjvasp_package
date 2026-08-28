"""Reusable DOI metadata, journal abbreviation, and PDF filename helpers.

Functions:
    normalize_doi: Strip DOI URL and ``doi:`` prefixes.
    parse_dois: Read a UTF-8 DOI list while preserving order.
    fetch_doi_metadata: Fetch naming metadata from Crossref.
    get_journal_abbreviation: Resolve a journal name through the shared index.
    check_journal_metadata: Reject non-journal and unsupported records.
    generate_pdf_filename: Build ``year-journal-title.pdf`` names.
    is_complete_pdf: Check PDF magic bytes and terminal EOF marker.
"""

from __future__ import annotations

import html
import json
import re
from pathlib import Path
from typing import Any
from urllib.parse import quote
from urllib.request import Request, urlopen


JOURNAL_ABBREVIATIONS = {
    "Acta Materialia": "ACTA-MATER",
    "Advanced Materials": "ADV-MATER",
    "Annual Review of Psychology": "ANNU-REV-PSYCHOL",
    "Applied Physics Letters": "APL",
    "BMC Bioinformatics": "BMC-BIOINFORMATICS",
    "Cell": "CELL",
    "Chemical Society Reviews": "CHEM-SOC-REV",
    "Communications in Mathematical Physics": "COMMUN-MATH-PHYS",
    "Computational Materials Science": "COMP-MATER-SCI",
    "Frontiers in Neuroscience": "FRONT-NEUROSCI",
    "IEEE Transactions on Information Theory": "IEEE-TIT",
    "Journal of the American Chemical Society": "JACS",
    "Journal of Physics: Condensed Matter": "JPCM",
    "Nature": "NATURE",
    "Nucleic Acids Research": "NAR",
    "Physical Review B": "PRB",
    "PLOS ONE": "PLOS-ONE",
    "Proceedings of the National Academy of Sciences": "PNAS",
    "Science": "SCIENCE",
    "Sensors": "SENSORS",
}

NON_ARTICLE_TITLE_PREFIXES = ("snapshot:",)


def normalize_doi(raw: str) -> str:
    """Return a DOI without URL or ``doi:`` prefixes."""
    value = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", raw.strip(), flags=re.I)
    return re.sub(r"^doi:\s*", "", value, flags=re.I).strip()


def parse_dois(path_file: Path) -> list[str]:
    """Read and de-duplicate the first field of each DOI-list line.

    Args:
        path_file: UTF-8 or UTF-8-with-BOM DOI list.

    Returns:
        DOI values in their first-seen order.
    """
    ldoi = []
    seen = set()
    for line in path_file.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        doi = normalize_doi(line.split()[0])
        if doi and doi not in seen:
            seen.add(doi)
            ldoi.append(doi)
    return ldoi


def fetch_doi_metadata(
        doi: str,
        timeout: int = 20,
        user_agent: str = "mymetal-literature/1.0") -> dict[str, Any] | None:
    """Fetch Crossref metadata used only for validation and naming.

    Args:
        doi: DOI value.
        timeout: HTTP timeout in seconds.
        user_agent: User-Agent sent to Crossref.

    Returns:
        Crossref's ``message`` object, or ``None`` when unavailable.
    """
    request = Request(
        "https://api.crossref.org/works/" + quote(normalize_doi(doi), safe=""),
        headers={"Accept": "application/json", "User-Agent": user_agent},
    )
    try:
        with urlopen(request, timeout=timeout) as response:
            return json.load(response).get("message")
    except (OSError, ValueError, json.JSONDecodeError):
        return None


def _get_first(value: Any) -> str:
    if isinstance(value, list):
        value = value[0] if value else ""
    return html.unescape(str(value or "")).strip()


def _normalize_journal_name(name: str) -> str:
    return "".join(char.lower() for char in name if char.isalnum())


def get_journal_abbreviation(
        journal: str,
        dict_journal_abbreviations: dict[str, str] | None = None) -> str | None:
    """Resolve a full or already abbreviated journal name.

    Args:
        journal: Journal name from publication metadata.
        dict_journal_abbreviations: Optional replacement/extension index.

    Returns:
        Uppercase filename-safe abbreviation, or ``None`` if unsupported.
    """
    dict_abbreviations = dict(JOURNAL_ABBREVIATIONS)
    if dict_journal_abbreviations:
        dict_abbreviations.update(dict_journal_abbreviations)
    normalized = _normalize_journal_name(journal)
    for name, abbreviation in dict_abbreviations.items():
        if normalized in {
                _normalize_journal_name(name),
                _normalize_journal_name(abbreviation)}:
            return re.sub(r"[^A-Za-z0-9]+", "-", abbreviation).strip("-").upper()
    return None


def get_publication_year(dict_metadata: dict[str, Any]) -> str:
    """Extract the best available four-digit publication year."""
    for key in ("published-print", "published-online", "issued", "created"):
        value = dict_metadata.get(key, {})
        if not isinstance(value, dict):
            continue
        lparts = value.get("date-parts", [[]])
        if lparts and lparts[0]:
            year = str(lparts[0][0])
            if re.fullmatch(r"\d{4}", year):
                return year
    year = str(dict_metadata.get("year", ""))
    return year if re.fullmatch(r"\d{4}", year) else ""


def check_journal_metadata(
        dict_metadata: dict[str, Any] | None,
        dict_journal_abbreviations: dict[str, str] | None = None) -> str | None:
    """Return why a record cannot be named as a supported journal article."""
    if not dict_metadata:
        return "metadata unavailable"
    if dict_metadata.get("type") != "journal-article":
        return "not a journal article"
    title = _get_first(dict_metadata.get("title"))
    if any(title.lower().startswith(prefix) for prefix in NON_ARTICLE_TITLE_PREFIXES):
        return "excluded article type"
    journal = _get_first(dict_metadata.get("container-title"))
    if not get_journal_abbreviation(journal, dict_journal_abbreviations):
        return "journal abbreviation missing"
    if not get_publication_year(dict_metadata):
        return "publication year missing"
    if not title:
        return "title missing"
    return None


def _get_title_prefix(title: str, character_count: int) -> str:
    title = re.sub(r"<[^>]+>", " ", html.unescape(title))
    lparts = []
    remaining = character_count
    for token in re.findall(r"[A-Za-z0-9]+|[\u3400-\u9fff]+", title):
        lparts.append(token[:remaining])
        remaining -= len(lparts[-1])
        if remaining <= 0:
            break
    return "-".join(lparts)


def generate_pdf_filename(
        dict_metadata: dict[str, Any],
        title_character_count: int = 10,
        dict_journal_abbreviations: dict[str, str] | None = None) -> str:
    """Build ``year-journal-first-ten-title-characters.pdf``.

    Args:
        dict_metadata: Crossref-style publication metadata.
        title_character_count: Maximum letters, digits, or Chinese characters.
        dict_journal_abbreviations: Optional replacement/extension index.

    Returns:
        Filename ending in ``.pdf``.

    Raises:
        ValueError: Metadata is unsupported or incomplete.
    """
    reason = check_journal_metadata(dict_metadata, dict_journal_abbreviations)
    if reason:
        raise ValueError(reason)
    title = _get_first(dict_metadata.get("title"))
    title_prefix = _get_title_prefix(title, title_character_count)
    if not title_prefix:
        raise ValueError("title has no filename-safe tokens")
    journal = _get_first(dict_metadata.get("container-title"))
    abbreviation = get_journal_abbreviation(journal, dict_journal_abbreviations)
    year = get_publication_year(dict_metadata)
    return "-".join([year, abbreviation, title_prefix]) + ".pdf"


def is_complete_pdf(path_pdf: Path, minimum_size: int = 5000) -> bool:
    """Return whether a PDF has a header, terminal EOF marker, and useful size."""
    try:
        size = path_pdf.stat().st_size
        if size < minimum_size:
            return False
        with path_pdf.open("rb") as file_pdf:
            if file_pdf.read(5) != b"%PDF-":
                return False
            file_pdf.seek(max(0, size - 2048))
            return b"%%EOF" in file_pdf.read()
    except OSError:
        return False
