"""FastAPI 入口。

两条相互独立的流程：
1. /compare：上传单条染色体参考 FASTA + 两份单样本 VCF，展示归一化核对结果；
2. /convert：上传旧、新单染色体 FASTA、旧参考单样本 VCF 与 TSV 坐标映射表，
   把按旧参考调用的变异整批转换到新参考，提供带源行号与 GT 的目标 VCF 下载。
"""

from pathlib import Path

from fastapi import FastAPI, File, Request, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse, Response
from fastapi.templating import Jinja2Templates

from .bio import (
    MAX_FASTA_BYTES,
    MAX_VCF_BYTES,
    CheckError,
    compare,
    parse_fasta,
    parse_vcf,
)
from .convert import MAX_MAP_BYTES, convert_vcf

app = FastAPI(
    title="VCF 变异工具集",
    description=(
        "双 VCF 归一化核对（/compare），以及旧参考到新参考的单样本 VCF 坐标转换（/convert）。"
    ),
)

templates = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))


def _read_limited(file: UploadFile, limit: int, label: str) -> bytes:
    data = file.file.read(limit + 1)
    if len(data) > limit:
        raise CheckError(
            f"{label} 超过大小上限 {limit} 字节"
            f"（FASTA 上限 {MAX_FASTA_BYTES // 1024} KiB，"
            f"VCF/映射表上限 {MAX_VCF_BYTES // 1024} KiB）。"
        )
    if not data:
        raise CheckError(f"{label} 为空文件。")
    return data


def _run_compare(fasta_bytes: bytes, vcf_a_bytes: bytes, vcf_b_bytes: bytes):
    refobj = parse_fasta(fasta_bytes)
    calls_a, sample_a = parse_vcf(vcf_a_bytes, "VCF A", refobj)
    calls_b, sample_b = parse_vcf(vcf_b_bytes, "VCF B", refobj)
    groups = compare(refobj, calls_a, calls_b)
    return refobj, sample_a, sample_b, groups


@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    return templates.TemplateResponse(
        request,
        "index.html",
        {
            "fasta_limit_kib": MAX_FASTA_BYTES // 1024,
            "vcf_limit_kib": MAX_VCF_BYTES // 1024,
        },
    )


@app.post("/compare", response_class=HTMLResponse)
def compare_page(
    request: Request,
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        refobj, sample_a, sample_b, groups = _run_compare(
            fasta_bytes, vcf_a_bytes, vcf_b_bytes
        )
    except CheckError as exc:
        # 整次核对失败：只展示错误，不展示任何局部结论
        return templates.TemplateResponse(
            request,
            "index.html",
            {
                "fasta_limit_kib": MAX_FASTA_BYTES // 1024,
                "vcf_limit_kib": MAX_VCF_BYTES // 1024,
                "error": str(exc),
            },
            status_code=422,
        )

    return templates.TemplateResponse(
        request,
        "result.html",
        {
            "ref_name": refobj.name,
            "ref_length": len(refobj.sequence),
            "sample_a": sample_a,
            "sample_b": sample_b,
            "groups": groups,
        },
    )


@app.post("/api/compare")
def compare_api(
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    """JSON 接口，便于脚本化调用；失败返回 422 与错误信息。"""
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        refobj, sample_a, sample_b, groups = _run_compare(
            fasta_bytes, vcf_a_bytes, vcf_b_bytes
        )
    except CheckError as exc:
        return JSONResponse(status_code=422, content={"ok": False, "error": str(exc)})

    def dump_entry(entry):
        def dump_call(call):
            if call is None:
                return None
            return {
                "dosage": call.dosage,
                "line_numbers": call.line_numbers,
                "raw": call.raw_descriptions,
            }

        v = entry["variant"]
        return {
            "variant": {
                "pos": v.pos,
                "ref": v.ref,
                "alt": v.alt,
                "kind": v.kind,
                "label": v.label,
            },
            "left": dump_call(entry["left"]),
            "right": dump_call(entry["right"]),
        }

    return {
        "ok": True,
        "reference": {"name": refobj.name, "length": len(refobj.sequence)},
        "samples": {"left": sample_a, "right": sample_b},
        "result": {key: [dump_entry(e) for e in rows] for key, rows in groups.items()},
    }


# --------------------------------------------------------------------------- #
# 旧参考 → 新参考 VCF 坐标转换
# --------------------------------------------------------------------------- #
def _index_context(**extra):
    ctx = {
        "fasta_limit_kib": MAX_FASTA_BYTES // 1024,
        "vcf_limit_kib": MAX_VCF_BYTES // 1024,
        "map_limit_kib": MAX_MAP_BYTES // 1024,
    }
    ctx.update(extra)
    return ctx


def _run_convert(old_fasta: bytes, new_fasta: bytes, mapping: bytes, vcf: bytes, name: str) -> bytes:
    return convert_vcf(old_fasta, new_fasta, mapping, vcf, vcf_filename=name)


def _vcf_record_count(target: bytes) -> int:
    return sum(1 for ln in target.splitlines() if ln and not ln.startswith(b"#"))


def _read_convert_inputs(
    old_fasta: UploadFile,
    new_fasta: UploadFile,
    mapping: UploadFile,
    vcf: UploadFile,
) -> bytes:
    old_bytes = _read_limited(old_fasta, MAX_FASTA_BYTES, "旧参考 FASTA")
    new_bytes = _read_limited(new_fasta, MAX_FASTA_BYTES, "新参考 FASTA")
    map_bytes = _read_limited(mapping, MAX_MAP_BYTES, "坐标映射表")
    vcf_bytes = _read_limited(vcf, MAX_VCF_BYTES, "源 VCF")
    return _run_convert(old_bytes, new_bytes, map_bytes, vcf_bytes, vcf.filename or "source.vcf")


def _vcf_response(target: bytes) -> Response:
    return Response(
        content=target,
        media_type="text/vcf; charset=utf-8",
        headers={
            "Content-Disposition": 'attachment; filename="target.vcf"',
            "X-Record-Count": str(_vcf_record_count(target)),
        },
    )


@app.get("/convert", response_class=HTMLResponse)
def convert_page(request: Request):
    return templates.TemplateResponse(request, "convert.html", _index_context())


@app.post("/convert", response_class=HTMLResponse)
def convert_submit(
    request: Request,
    old_fasta: UploadFile = File(...),
    new_fasta: UploadFile = File(...),
    mapping: UploadFile = File(...),
    vcf: UploadFile = File(...),
):
    """浏览器表单入口：成功直接下载目标 VCF，失败返回整批错误页（无部分文件）。"""
    try:
        target = _read_convert_inputs(old_fasta, new_fasta, mapping, vcf)
    except CheckError as exc:
        return templates.TemplateResponse(
            request,
            "convert.html",
            _index_context(error=str(exc)),
            status_code=422,
        )
    return _vcf_response(target)


@app.post("/api/convert")
def convert_api(
    old_fasta: UploadFile = File(...),
    new_fasta: UploadFile = File(...),
    mapping: UploadFile = File(...),
    vcf: UploadFile = File(...),
):
    """JSON 接口：成功直接返回目标 VCF 文件流，失败返回 422 与错误信息。"""
    try:
        target = _read_convert_inputs(old_fasta, new_fasta, mapping, vcf)
    except CheckError as exc:
        return JSONResponse(status_code=422, content={"ok": False, "error": str(exc)})
    return _vcf_response(target)
