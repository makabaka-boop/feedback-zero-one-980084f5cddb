"""非法输入验收：整次拒绝、上报偏移正确、输出文件保持原样。"""

import os

import pytest

from conftest import run_cli


def error_entries(manifest):
    return sorted(k for k, v in manifest.items() if v["status"] == "error")


@pytest.mark.parametrize("name", [
    # BGZF 层：报压缩文件零基块起点
    "bgzf_empty_file", "bgzf_bad_magic", "bgzf_bad_cm", "bgzf_flg0",
    "bgzf_flg6", "bgzf_two_bc", "bgzf_no_bc", "bgzf_bc_slen4",
    "bgzf_bad_crc", "bgzf_bad_isize", "bgzf_bad_deflate",
    "bgzf_bsize_small", "bgzf_bsize_overrun", "bgzf_isize_over_64k",
    "bgzf_truncated_tail", "bgzf_no_eof", "bgzf_garbage_after_eof",
    "bgzf_second_block_crc",
    # 头部：报 0
    "hdr_bad_magic", "hdr_short", "hdr_ltext_negative", "hdr_ltext_overrun",
    "hdr_nref_negative", "hdr_nref_overrun", "hdr_ref_dup",
    "hdr_ref_len_zero", "hdr_ref_len_negative", "hdr_ref_name_no_nul",
    "hdr_ref_lname_zero", "hdr_ref_truncated",
    # 记录：报 vo
    "rec_bad_op", "rec_op_b", "rec_zero_len", "rec_consumption", "rec_span",
    "rec_pos_beyond", "rec_soft_middle", "rec_hard_middle",
    "rec_h_before_s_end", "rec_qual_mixed", "rec_qual_94", "rec_dup_aux",
    "rec_bad_aux_type", "rec_aux_truncated", "rec_aux_b_bad_sub",
    "rec_aux_b_overrun", "rec_aux_extra_byte", "rec_placeholder_no_cg",
    "rec_placeholder_cg_z", "rec_refid_range", "rec_pos_negative",
    "rec_name_no_nul", "rec_block_size_small", "rec_block_size_overrun",
    "rec_leftover", "rec_second_record_bad_multiblock",
])
def test_rejected_with_offset(name, manifest, fixtures, tmp_path):
    entry = manifest[name]
    assert entry["status"] == "error"
    out = tmp_path / "out.json"
    out.write_text("OLD")  # 预置旧文件：失败时必须保留
    res = run_cli(os.path.join(fixtures, entry["bam"]), out)
    assert res.returncode == 1, (name, res.returncode, res.stdout, res.stderr)
    assert res.stdout == ""
    assert res.stderr.strip() == str(entry["offset"]), (
        name, res.stderr, entry["offset"])
    assert out.read_text() == "OLD"
    assert not list(tmp_path.glob("*.tmp.*"))  # 无残留临时文件


def test_error_offsets_are_plausible(manifest):
    """偏移分类自检：头错恒 0；跨块记录错的 vo 含非零压缩起点。"""
    for k, v in manifest.items():
        if v["status"] != "error":
            continue
        if k.startswith("hdr_"):
            assert v["offset"] == 0
    assert manifest["rec_second_record_bad_multiblock"]["offset"] > 0x10000
    assert manifest["bgzf_second_block_crc"]["offset"] > 0
    assert manifest["bgzf_no_eof"]["offset"] > 0
