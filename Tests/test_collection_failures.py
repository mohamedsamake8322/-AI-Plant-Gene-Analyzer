import requests
import pytest

from collect import collect_kegg
import request_utils
from scripts import collect_multi_type


def test_http_helper_does_not_start_a_second_unbounded_retry_loop(monkeypatch):
    calls = {"session": 0, "fallback": 0}

    class FailingSession:
        def get(self, *args, **kwargs):
            calls["session"] += 1
            raise requests.Timeout("probe timeout")

    def fallback_get(*args, **kwargs):
        calls["fallback"] += 1
        raise AssertionError("raw requests fallback must not be used")

    monkeypatch.setattr(request_utils, "_throttle", lambda _url: None)
    monkeypatch.setattr(request_utils, "get_session", lambda retries=2: FailingSession())
    monkeypatch.setattr(request_utils.requests, "get", fallback_get)

    with pytest.raises(requests.Timeout, match="probe timeout"):
        request_utils.get("https://example.invalid/resource", timeout=1)

    assert calls == {"session": 1, "fallback": 0}


def test_kegg_gene_list_transport_failure_is_not_reported_as_empty_data(monkeypatch):
    def fail_request(*args, **kwargs):
        raise requests.Timeout("KEGG unavailable")

    monkeypatch.setattr(collect_kegg.rq, "get", fail_request)

    with pytest.raises(RuntimeError, match="KEGG gene-list request failed"):
        collect_kegg._get_gene_list("zma", 5)


def test_ncbi_sequence_type_failure_can_be_reported_to_orchestrator(monkeypatch, tmp_path):
    def fail_pipeline(_argv):
        raise requests.Timeout("NCBI timeout")

    monkeypatch.setattr(collect_multi_type.pipeline_module, "main", fail_pipeline)

    with pytest.raises(RuntimeError, match="NCBI protein collection failed"):
        collect_multi_type.collect_and_clean_type(
            "Zea mays",
            "protein",
            5,
            tmp_path / "raw.json",
            tmp_path / "clean.json",
            raise_on_error=True,
        )
