from app.services import omim_store, panel_deadzone


def _row(symbol, disease):
    return {"gene_symbol": symbol, "OMIM_id": "", "OMIM_disease": disease}


def test_disease_associated_includes_platform_omim_disease_genes(monkeypatch):
    rows = {
        "TBPL2": _row("TBPL2", "Oocyte/zygote/embryo maturation arrest 29 (621770)(AR)"),
        "SUSCX1": _row("SUSCX1", "{Susceptibility to something} (622001)(AD)"),
        "PROVX1": _row("PROVX1", "?Provisional disorder (622002)(AR)"),
        "BLOODX1": _row("BLOODX1", "[Blood group, made up] (622003)"),
        "NOMIMX1": _row("NOMIMX1", "Disease without MIM (AR)"),
    }
    monkeypatch.setattr(omim_store, "_state", {"mtime": ("test", 1), "by_gene": rows, "by_omim_id": {}})
    monkeypatch.setitem(panel_deadzone._omim_cache, "key", None)
    assert panel_deadzone.is_disease_associated_gene("TBPL2")
    assert panel_deadzone.is_disease_associated_gene("SUSCX1")
    assert panel_deadzone.is_disease_associated_gene("PROVX1")
    assert not panel_deadzone.is_disease_associated_gene("BLOODX1")
    assert not panel_deadzone.is_disease_associated_gene("NOMIMX1")
    # Delivered panel genes still count; dead-zone tables only cover them.
    assert panel_deadzone.is_disease_associated_gene("BRCA1")
    assert panel_deadzone.has_dead_zone_evaluation("BRCA1")
    assert not panel_deadzone.has_dead_zone_evaluation("TBPL2")
    # A new OMIM revision is picked up without restarting.
    rows2 = {"TBPL2": _row("TBPL2", "Name only")}
    monkeypatch.setattr(omim_store, "_state", {"mtime": ("test", 2), "by_gene": rows2, "by_omim_id": {}})
    assert not panel_deadzone.is_disease_associated_gene("TBPL2")
    monkeypatch.setitem(panel_deadzone._omim_cache, "key", None)
