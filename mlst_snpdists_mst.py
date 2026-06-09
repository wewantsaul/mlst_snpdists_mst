#!/usr/bin/env python3
"""
mlst_snpdists_mst.py
====================
Generate an interactive Minimum Spanning Tree (MST) HTML visualizer
from MLST (tseemann) and snp-dists outputs.

by: Saul M. Rojas

Usage
-----
  python mlst_snpdists_mst.py --mlst mlst.txt --snpdists snp-dists.txt
  python mlst_snpdists_mst.py --mlst mlst.txt --snpdists snp-dists.txt --output my_tree.html
  python mlst_snpdists_mst.py --mlst mlst.txt --snpdists snp-dists.txt --title "E. coli outbreak 2025"

Supports
--------
  - Legacy MLST format  (FILE  SCHEME  ST  adk  fumC ...)
  - Non-legacy format   (sample  scheme  ST  adk(53)  fumC(40) ...)
  Format is auto-detected.

  Samples with no ST ("-") are shown as grey nodes.

Output
------
  A single self-contained HTML file with no external dependencies.
  The page lets you:
    - Adjust SNP threshold to define clusters (MST redraws live)
    - Pan (drag) and zoom (scroll wheel)
    - Hover nodes for allele profile tooltips
    - Download the current view as PNG or JPEG
"""

import argparse
import json
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def _norm_name(s: str) -> str:
    """Strip .fasta / .fasta.ref suffixes used by parsnp / snp-dists."""
    return re.sub(r"\.fasta(\.ref)?$", "", s.strip(), flags=re.IGNORECASE)


def parse_mlst(path: str) -> dict:
    """
    Parse MLST output (legacy or no-legacy).

    Returns
    -------
    dict  {sample_name: {"st": str, "scheme": str, "alleles": {gene: allele}}}
    """
    lines = Path(path).read_text().strip().splitlines()
    lines = [l for l in lines if l.strip()]
    if not lines:
        sys.exit(f"ERROR: MLST file '{path}' is empty.")

    result = {}
    first_cols = lines[0].split("\t")

    # Legacy format: first column header is FILE or SAMPLE (case-insensitive)
    is_legacy = first_cols[0].strip().upper() in ("FILE", "SAMPLE")

    if is_legacy:
        header = first_cols
        gene_start = 3
        genes = [h.strip() for h in header[gene_start:]]
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            name = _norm_name(parts[0])
            st = parts[2].strip() if parts[2].strip() not in ("?", "") else "-"
            alleles = {}
            for j, gene in enumerate(genes):
                idx = gene_start + j
                alleles[gene] = parts[idx].strip() if idx < len(parts) else "?"
            result[name] = {"st": st, "scheme": parts[1].strip(), "alleles": alleles}
    else:
        # No-legacy format: gene(allele) columns, no header row
        for line in lines:
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            name = _norm_name(parts[0])
            st = parts[2].strip() if parts[2].strip() not in ("?", "") else "-"
            alleles = {}
            for col in parts[3:]:
                m = re.match(r"^(\w+)\((.+)\)$", col.strip())
                if m:
                    alleles[m.group(1)] = m.group(2)
            result[name] = {"st": st, "scheme": parts[1].strip(), "alleles": alleles}

    return result


def parse_snpdists(path: str) -> dict:
    """
    Parse snp-dists tab-separated matrix.

    Returns
    -------
    {"names": [str, ...], "matrix": {row: {col: int}}}
    """
    lines = Path(path).read_text().strip().splitlines()
    lines = [l for l in lines if l.strip()]
    if not lines:
        sys.exit(f"ERROR: snp-dists file '{path}' is empty.")

    # snp-dists first line: "snp-dists 1.x.x\tSample1\tSample2\t..."
    # The version token and column names are on the SAME line, tab-separated.
    # Data rows start from line index 1.
    header = lines[0].split("\t")
    # If first token looks like a version string, column names start at index 1
    if re.match(r"snp-dists", header[0], re.IGNORECASE):
        col_names = [_norm_name(c) for c in header[1:]]
        data_start = 1
    else:
        # Plain matrix with no version header
        col_names = [_norm_name(c) for c in header[1:]]
        data_start = 1

    matrix = {}
    for line in lines[data_start:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        row_name = _norm_name(parts[0])
        matrix[row_name] = {}
        for j, col in enumerate(col_names):
            try:
                matrix[row_name][col] = int(parts[j + 1])
            except (ValueError, IndexError):
                matrix[row_name][col] = 0

    # De-duplicate column names (e.g. EC02.fasta and EC02.fasta.ref collapse)
    seen = set()
    unique_names = []
    for n in col_names:
        if n not in seen:
            seen.add(n)
            unique_names.append(n)

    return {"names": unique_names, "matrix": matrix}


def get_max_snp(snp_data: dict) -> int:
    mx = 0
    for row in snp_data["matrix"].values():
        for v in row.values():
            if v > mx:
                mx = v
    return mx


# ---------------------------------------------------------------------------
# HTML template (self-contained, no external deps)
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{page_title}</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Arial,sans-serif;
     background:#f5f5f0;color:#1a1a18;font-size:14px}}
header{{background:#fff;border-bottom:1px solid #e0dfd7;padding:12px 20px;
        display:flex;align-items:center;gap:12px}}
header h1{{font-size:17px;font-weight:600;color:#1a1a18}}
header span{{font-size:12px;color:#6b6b65;border:1px solid #d0cfc7;
             border-radius:6px;padding:2px 8px}}
.main{{display:flex;height:calc(100vh - 49px)}}
.sidebar{{width:260px;flex-shrink:0;background:#fff;border-right:1px solid #e0dfd7;
          padding:16px;overflow-y:auto;display:flex;flex-direction:column;gap:16px}}
.sidebar h2{{font-size:13px;font-weight:600;color:#4a4a46;text-transform:uppercase;
             letter-spacing:.04em;margin-bottom:4px}}
.ctrl-group{{display:flex;flex-direction:column;gap:10px}}
label.ctrl-label{{font-size:12px;color:#6b6b65;display:flex;
                  justify-content:space-between;align-items:center}}
label.ctrl-label span{{font-weight:600;color:#1a1a18;font-size:13px}}
input[type=range]{{width:100%;accent-color:#378ADD;cursor:pointer}}
.stat-grid{{display:grid;grid-template-columns:1fr 1fr;gap:8px}}
.stat-card{{background:#f5f5f0;border-radius:8px;padding:8px 12px}}
.stat-card .sv{{font-size:22px;font-weight:600;line-height:1}}
.stat-card .sl{{font-size:11px;color:#6b6b65;margin-top:2px}}
.btn{{width:100%;padding:8px;cursor:pointer;border-radius:8px;
      border:1px solid #d0cfc7;background:#fff;color:#1a1a18;
      font-size:13px;text-align:left;display:flex;align-items:center;gap:6px}}
.btn:hover{{background:#f5f5f0}}
.btn svg{{flex-shrink:0}}
.legend-wrap{{display:flex;flex-direction:column;gap:6px}}
.leg-item{{display:flex;align-items:center;gap:7px;font-size:12px}}
.leg-dot{{width:12px;height:12px;border-radius:50%;flex-shrink:0}}
.warn{{background:#FAEEDA;border:1px solid #F0C98A;border-radius:8px;
       padding:8px 10px;font-size:12px;color:#854F0B;line-height:1.5}}
.canvas-wrap{{flex:1;position:relative;overflow:hidden;background:#fafaf7}}
canvas{{display:block;cursor:grab}}
canvas.panning{{cursor:grabbing}}
.tooltip{{position:absolute;background:#fff;border:1px solid #d0cfc7;
          border-radius:8px;padding:8px 12px;font-size:12px;pointer-events:none;
          opacity:0;transition:opacity .15s;max-width:230px;z-index:10;
          line-height:1.65;box-shadow:0 2px 8px rgba(0,0,0,.08)}}
.tooltip b{{font-weight:600;display:block;margin-bottom:3px;font-size:13px}}
.tooltip .al{{color:#6b6b65;font-size:11px;margin-top:4px;border-top:1px solid #e8e8e2;
              padding-top:4px}}
.footer{{font-size:11px;color:#9b9b96;padding:6px 12px;
         background:#fff;border-top:1px solid #e8e8e2;text-align:center}}
</style>
</head>
<body>
<header>
  <h1>{page_title}</h1>
  <span>MLST + SNP-dists MST</span>
</header>
<div class="main">
  <div class="sidebar">
    <div>
      <h2>Clustering</h2>
      <div class="ctrl-group" style="margin-top:8px">
        <label class="ctrl-label">SNP threshold</label>
        <div style="display:flex;align-items:center;gap:8px">
          <input type="range" id="thresh" min="0" max="{max_snp}" value="{default_thresh}" step="1" style="flex:1;min-width:0">
          <input type="number" id="thresh-num" min="0" max="{max_snp}" value="{default_thresh}"
                 style="width:72px;padding:4px 6px;border:1px solid #d0cfc7;border-radius:6px;
                        font-size:13px;font-family:inherit;text-align:right;flex-shrink:0">
        </div>
        <div style="font-size:11px;color:#9b9b96">Samples ≤ threshold SNPs apart are in the same cluster. Edges above threshold are hidden.</div>
      </div>
    </div>

    <div>
      <h2>Display</h2>
      <div class="ctrl-group" style="margin-top:8px">
        <label class="ctrl-label">Node size <span id="node-sz-val">22</span></label>
        <input type="range" id="node-sz" min="10" max="42" value="22" step="1">
        <label class="ctrl-label">Label size <span id="font-sz-val">11</span>px</label>
        <input type="range" id="font-sz" min="7" max="18" value="11" step="1">
        <label class="ctrl-label" style="gap:8px">
          <input type="checkbox" id="show-edge-labels" checked>
          Show SNP distances on edges
        </label>
      </div>
    </div>

    <div>
      <h2>Statistics</h2>
      <div class="stat-grid" id="stats-grid">
        <div class="stat-card"><div class="sv" id="s-samples">-</div><div class="sl">samples</div></div>
        <div class="stat-card"><div class="sv" id="s-sts">-</div><div class="sl">unique STs</div></div>
        <div class="stat-card"><div class="sv" id="s-clusters">-</div><div class="sl">clusters</div></div>
        <div class="stat-card"><div class="sv" id="s-no-st">-</div><div class="sl">no ST</div></div>
      </div>
    </div>

    <div id="warn-box" class="warn" style="display:none"></div>

    <div>
      <h2>Export</h2>
      <div style="display:flex;flex-direction:column;gap:6px;margin-top:8px">
        <button class="btn" onclick="saveImg('png')">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
          Save as PNG
        </button>
        <button class="btn" onclick="saveImg('jpeg')">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
          Save as JPEG
        </button>
        <button class="btn" onclick="resetView()">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><polyline points="1 4 1 10 7 10"/><path d="M3.51 15a9 9 0 1 0 .49-3.5"/></svg>
          Reset view
        </button>
      </div>
    </div>

    <div>
      <h2>Legend</h2>
      <div class="legend-wrap" id="legend" style="margin-top:8px"></div>
    </div>
  </div>

  <div class="canvas-wrap">
    <canvas id="cvs"></canvas>
    <div class="tooltip" id="tip"></div>
    <div style="position:absolute;bottom:8px;right:12px;font-size:11px;color:#9b9b96">
      Drag nodes to reposition &nbsp;·&nbsp; Drag background to pan &nbsp;·&nbsp; Scroll to zoom &nbsp;·&nbsp; Hover for details
    </div>
  </div>
</div>

<script>
const MLST_DATA = {mlst_json};
const SNP_NAMES = {snp_names_json};
const SNP_MATRIX = {snp_matrix_json};
const MAX_SNP = {max_snp};

const ST_COLORS = [
  '#378ADD','#1D9E75','#D85A30','#7F77DD','#BA7517',
  '#D4537E','#639922','#E24B4A','#0F6E56','#533AB7',
  '#993C1D','#188FA7','#C06014','#8A4FFF','#2D7D46'
];
const NO_ST_COLOR = '#888780';

let thresh = {default_thresh};
let nodeR = 22, fontSize = 11, showEdgeLabels = true;
let pan = {{x:0,y:0}}, scale = 1;
let dragging = false, dragStart = {{x:0,y:0}}, panStart = {{x:0,y:0}};

let nodes = [], mstEdges = [], clusters = [], nodePos = {{}};
let stColorMap = {{}};
let cvs, ctx, tip, W, H;

function normName(s){{ return s.replace(/\.fasta(\.ref)?$/i,'').trim(); }}

function buildMST(names, matrix){{
  const edgeList = [];
  for(let i=0;i<names.length;i++) for(let j=i+1;j<names.length;j++){{
    const a=names[i], b=names[j];
    let d = Infinity;
    if(matrix[a]&&matrix[a][b]!=null) d=matrix[a][b];
    else if(matrix[b]&&matrix[b][a]!=null) d=matrix[b][a];
    edgeList.push({{a,b,d}});
  }}
  edgeList.sort((x,y)=>x.d-y.d);
  const parent={{}};
  names.forEach(n=>parent[n]=n);
  function find(x){{return parent[x]===x?x:parent[x]=find(parent[x]);}}
  function union(x,y){{parent[find(x)]=find(y);}}
  const mst=[];
  for(const e of edgeList){{
    if(find(e.a)!==find(e.b)){{union(e.a,e.b);mst.push(e);}}
  }}
  return mst;
}}

function getClusters(names, mst, threshold){{
  const parent={{}};
  names.forEach(n=>parent[n]=n);
  function find(x){{return parent[x]===x?x:parent[x]=find(parent[x]);}}
  function union(x,y){{parent[find(x)]=find(y);}}
  for(const e of mst){{if(e.d<=threshold) union(e.a,e.b);}}
  const groups={{}};
  names.forEach(n=>{{const r=find(n);if(!groups[r])groups[r]=[];groups[r].push(n);}});
  return Object.values(groups);
}}

function layoutForce(names, mst, pos, lW, lH, pad){{
  const n = names.length;
  const r = (Math.min(lW,lH)/2 - pad) * 0.88;
  names.forEach((nm,i)=>{{
    const a = 2*Math.PI*i/n - Math.PI/2;
    pos[nm] = {{x: lW/2+r*Math.cos(a), y: lH/2+r*Math.sin(a)}};
  }});
  if(n<=1) return;
  const ideal = Math.min(lW,lH) / Math.max(n,2) * 1.4;
  for(let iter=0;iter<350;iter++){{
    const force={{}};
    names.forEach(nm=>force[nm]={{x:0,y:0}});
    for(let i=0;i<names.length;i++) for(let j=i+1;j<names.length;j++){{
      const a=names[i],b=names[j];
      const dx=pos[b].x-pos[a].x, dy=pos[b].y-pos[a].y;
      const d=Math.sqrt(dx*dx+dy*dy)||1;
      const rep=ideal*ideal/d;
      force[a].x-=rep*dx/d; force[a].y-=rep*dy/d;
      force[b].x+=rep*dx/d; force[b].y+=rep*dy/d;
    }}
    mst.forEach(e=>{{
      const dx=pos[e.b].x-pos[e.a].x, dy=pos[e.b].y-pos[e.a].y;
      const d=Math.sqrt(dx*dx+dy*dy)||1;
      const att=(d-ideal)*0.06;
      force[e.a].x+=att*dx/d; force[e.a].y+=att*dy/d;
      force[e.b].x-=att*dx/d; force[e.b].y-=att*dy/d;
    }});
    const cool=1-iter/350;
    names.forEach(nm=>{{
      pos[nm].x=Math.max(pad+nodeR, Math.min(lW-pad-nodeR, pos[nm].x+force[nm].x*cool*0.45));
      pos[nm].y=Math.max(pad+nodeR, Math.min(lH-pad-nodeR, pos[nm].y+force[nm].y*cool*0.45));
    }});
  }}
}}

function layoutMST(names, mst, threshold){{
  const pos={{}};
  clusters = getClusters(names, mst, threshold);
  if(clusters.length===1){{
    layoutForce(names, mst, pos, W, H, 60);
    return pos;
  }}
  const cols = Math.ceil(Math.sqrt(clusters.length));
  const rows = Math.ceil(clusters.length/cols);
  const cellW=(W-40)/cols, cellH=(H-40)/rows;
  clusters.forEach((cl,ci)=>{{
    const cx=20+cellW*(ci%cols)+cellW/2;
    const cy=20+cellH*Math.floor(ci/cols)+cellH/2;
    if(cl.length===1){{pos[cl[0]]={{x:cx,y:cy}};return;}}
    const subMst=mst.filter(e=>cl.includes(e.a)&&cl.includes(e.b));
    const subPos={{}};
    layoutForce(cl, subMst, subPos, cellW*0.88, cellH*0.88, 20);
    const cx0=Object.values(subPos).reduce((s,p)=>s+p.x,0)/cl.length;
    const cy0=Object.values(subPos).reduce((s,p)=>s+p.y,0)/cl.length;
    cl.forEach(n=>{{pos[n]={{x:cx+(subPos[n].x-cx0), y:cy+(subPos[n].y-cy0)}};}} );
  }});
  return pos;
}}

function draw(offCtx, offW, offH, bgColor){{
  const c = offCtx||ctx;
  const w = offW||W, h = offH||H;
  c.clearRect(0,0,w,h);
  if(bgColor){{c.fillStyle=bgColor;c.fillRect(0,0,w,h);}}
  c.save();
  const sc = offCtx ? scale : scale;
  const px = offCtx ? pan.x*(offW/W) : pan.x;
  const py = offCtx ? pan.y*(offH/H) : pan.y;
  c.translate(px,py);
  c.scale(sc,sc);

  const visEdges = mstEdges.filter(e=>e.d<=thresh);
  visEdges.forEach(e=>{{
    const a=nodePos[e.a], b=nodePos[e.b];
    if(!a||!b) return;
    c.beginPath();c.moveTo(a.x,a.y);c.lineTo(b.x,b.y);
    c.strokeStyle='rgba(136,135,128,0.4)';
    c.lineWidth=1.8/sc;c.stroke();
    if(showEdgeLabels){{
      const mx=(a.x+b.x)/2, my=(a.y+b.y)/2;
      c.fillStyle='rgba(90,90,85,0.85)';
      c.font=`${{Math.max(8,fontSize-2)/sc}}px -apple-system,Arial,sans-serif`;
      c.textAlign='center';c.textBaseline='middle';
      c.fillText(e.d.toLocaleString(), mx, my-8/sc);
    }}
  }});

  nodes.forEach(nd=>{{
    const p=nodePos[nd.id]; if(!p) return;
    const r=nodeR/sc;
    const col=stColorMap[nd.st]||NO_ST_COLOR;
    c.beginPath();c.arc(p.x,p.y,r,0,2*Math.PI);
    c.fillStyle=col;c.fill();
    c.strokeStyle='rgba(255,255,255,0.75)';c.lineWidth=2/sc;c.stroke();
    c.fillStyle='#fff';
    c.font=`${{fontSize/sc}}px -apple-system,Arial,sans-serif`;
    c.textAlign='center';c.textBaseline='middle';
    c.fillText(nd.id, p.x, p.y);
    const stLabel='ST'+(nd.st==='-'?'?':nd.st);
    c.fillStyle='rgba(255,255,255,0.9)';
    c.font=`${{Math.max(8,fontSize-1)/sc}}px -apple-system,Arial,sans-serif`;
    c.fillText(stLabel, p.x, p.y+r+10/sc);
  }});
  c.restore();
}}

function updateStats(){{
  const noST=nodes.filter(n=>n.st==='-').length;
  const stCount=new Set(nodes.filter(n=>n.st!=='-').map(n=>n.st)).size;
  document.getElementById('s-samples').textContent=nodes.length;
  document.getElementById('s-sts').textContent=stCount;
  document.getElementById('s-clusters').textContent=clusters.length;
  document.getElementById('s-no-st').textContent=noST;
  const wb=document.getElementById('warn-box');
  if(noST>0){{
    wb.style.display='';
    wb.textContent=noST+' sample'+(noST>1?'s':'')+' lack an ST assignment (shown in grey). '
      +'This typically indicates a novel allele or insufficient assembly quality — '
      +'they are still placed in the MST by SNP distance.';
  }} else wb.style.display='none';
}}

function buildLegend(){{
  const leg=document.getElementById('legend');
  leg.innerHTML='';
  const stList=[...new Set(nodes.map(n=>n.st))].sort((a,b)=>{{
    if(a==='-') return 1; if(b==='-') return -1;
    return parseInt(a)-parseInt(b)||a.localeCompare(b);
  }});
  stList.forEach(st=>{{
    const col=stColorMap[st]||NO_ST_COLOR;
    const div=document.createElement('div');
    div.className='leg-item';
    div.innerHTML=`<div class="leg-dot" style="background:${{col}}"></div><span>${{st==='-'?'No ST assigned':'ST '+st}}</span>`;
    leg.appendChild(div);
  }});
}}

function getNodeAt(mx,my){{
  const tx=(mx-pan.x)/scale, ty=(my-pan.y)/scale;
  for(const nd of nodes){{
    const p=nodePos[nd.id]; if(!p) continue;
    if(Math.hypot(tx-p.x,ty-p.y)<nodeR/scale) return nd;
  }}
  return null;
}}

function resetView(){{pan={{x:0,y:0}};scale=1;draw();}}

function saveImg(fmt){{
  const ratio=2;
  const oc=document.createElement('canvas');
  oc.width=W*ratio;oc.height=H*ratio;
  const ox=oc.getContext('2d');
  ox.scale(ratio,ratio);
  draw(ox,W*ratio,H*ratio,'#ffffff');
  const mime=fmt==='jpeg'?'image/jpeg':'image/png';
  oc.toBlob(blob=>{{
    const a=document.createElement('a');
    a.href=URL.createObjectURL(blob);
    a.download='mst_tree.'+fmt;
    a.click();
  }},mime,0.95);
}}

function init(){{
  cvs=document.getElementById('cvs');
  tip=document.getElementById('tip');
  ctx=cvs.getContext('2d');
  const wrap=cvs.parentElement;
  W=wrap.clientWidth;H=wrap.clientHeight;
  cvs.width=W;cvs.height=H;
  cvs.style.width=W+'px';cvs.style.height=H+'px';

  // Build node list from SNP names (source of truth for sample set)
  const stList=[...new Set(
    SNP_NAMES.map(n=>MLST_DATA[n]?MLST_DATA[n].st:'-')
  )].sort((a,b)=>{{if(a==='-')return 1;if(b==='-')return -1;return parseInt(a)-parseInt(b)||a.localeCompare(b);}});
  stList.forEach((st,i)=>{{stColorMap[st]=st==='-'?NO_ST_COLOR:ST_COLORS[i%ST_COLORS.length];}});

  nodes=SNP_NAMES.map(name=>{{
    const m=MLST_DATA[name]||{{}};
    return {{id:name, st:m.st||'-', scheme:m.scheme||'?', alleles:m.alleles||{{}}}};
  }});

  mstEdges=buildMST(SNP_NAMES,SNP_MATRIX);
  thresh=parseInt(document.getElementById('thresh').value);
  clusters=getClusters(SNP_NAMES,mstEdges,thresh);
  nodePos=layoutMST(SNP_NAMES,mstEdges,thresh);
  updateStats();
  buildLegend();
  draw();

  // Controls — threshold: slider and number input stay in sync
  function applyThresh(val){{
    thresh=Math.max(0,Math.min({max_snp},isNaN(val)?0:val));
    document.getElementById('thresh').value=thresh;
    document.getElementById('thresh-num').value=thresh;
    clusters=getClusters(SNP_NAMES,mstEdges,thresh);
    nodePos=layoutMST(SNP_NAMES,mstEdges,thresh);
    updateStats();draw();
  }}
  document.getElementById('thresh').addEventListener('input',function(){{
    applyThresh(parseInt(this.value));
  }});
  document.getElementById('thresh-num').addEventListener('input',function(){{
    applyThresh(parseInt(this.value));
  }});
  document.getElementById('thresh-num').addEventListener('keydown',function(e){{
    if(e.key==='ArrowUp'){{e.preventDefault();applyThresh(thresh+1);}}
    if(e.key==='ArrowDown'){{e.preventDefault();applyThresh(Math.max(0,thresh-1));}}
  }});
  document.getElementById('node-sz').addEventListener('input',function(){{
    nodeR=parseInt(this.value);
    document.getElementById('node-sz-val').textContent=nodeR;
    draw();
  }});
  document.getElementById('font-sz').addEventListener('input',function(){{
    fontSize=parseInt(this.value);
    document.getElementById('font-sz-val').textContent=fontSize;
    draw();
  }});
  document.getElementById('show-edge-labels').addEventListener('change',function(){{
    showEdgeLabels=this.checked;draw();
  }});

  // Interaction: drag node OR pan canvas
  let dragNode=null, isPanning=false;
  function canvasXY(e){{
    const r=cvs.getBoundingClientRect();
    return {{mx:e.clientX-r.left, my:e.clientY-r.top}};
  }}
  function toWorld(mx,my){{
    return {{wx:(mx-pan.x)/scale, wy:(my-pan.y)/scale}};
  }}
  cvs.addEventListener('mousedown',e=>{{
    const {{mx,my}}=canvasXY(e);
    const nd=getNodeAt(mx,my);
    if(nd){{
      dragNode=nd;
      tip.style.opacity='0';
      cvs.style.cursor='grabbing';
    }} else {{
      isPanning=true;
      cvs.classList.add('panning');
      dragStart={{x:mx,y:my}};
      panStart={{x:pan.x,y:pan.y}};
    }}
  }});
  window.addEventListener('mouseup',()=>{{
    dragNode=null;
    isPanning=false;
    cvs.classList.remove('panning');
    cvs.style.cursor='';
  }});
  cvs.addEventListener('mousemove',e=>{{
    const {{mx,my}}=canvasXY(e);
    if(dragNode){{
      const {{wx,wy}}=toWorld(mx,my);
      nodePos[dragNode.id]={{x:wx,y:wy}};
      draw();return;
    }}
    if(isPanning){{
      pan.x=panStart.x+(mx-dragStart.x);
      pan.y=panStart.y+(my-dragStart.y);
      draw();return;
    }}
    const nd=getNodeAt(mx,my);
    if(nd){{
      cvs.style.cursor='grab';
      tip.style.opacity='1';
      const alStr=Object.entries(nd.alleles)
        .map(([g,v])=>`${{g}}: ${{v}}`).join('  ·  ');
      tip.innerHTML=`<b>${{nd.id}}</b>ST: ${{nd.st==='–'||nd.st==='-'?'Not assigned':nd.st}}<br>Scheme: ${{nd.scheme}}`
        +(alStr?`<div class="al">${{alStr}}</div>`:'');
      tip.style.left=Math.min(mx+16,W-240)+'px';
      tip.style.top=Math.max(my-20,5)+'px';
    }} else {{
      cvs.style.cursor='';
      tip.style.opacity='0';
    }}
  }});
  cvs.addEventListener('mouseleave',()=>{{
    tip.style.opacity='0';
    if(!dragNode) cvs.style.cursor='';
  }});
  cvs.addEventListener('wheel',e=>{{
    e.preventDefault();
    const {{mx,my}}=canvasXY(e);
    const delta=e.deltaY>0?0.88:1.14;
    const ns=Math.min(6,Math.max(0.15,scale*delta));
    pan.x=mx-(mx-pan.x)*(ns/scale);
    pan.y=my-(my-pan.y)*(ns/scale);
    scale=ns;draw();
  }},{{passive:false}});

  window.addEventListener('resize',()=>{{
    W=wrap.clientWidth;H=wrap.clientHeight;
    cvs.width=W;cvs.height=H;
    cvs.style.width=W+'px';cvs.style.height=H+'px';
    nodePos=layoutMST(SNP_NAMES,mstEdges,thresh);
    draw();
  }});
}}

window.addEventListener('load',init);
</script>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate an interactive MST HTML visualizer from MLST + snp-dists outputs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("--mlst",      required=True, help="Path to MLST output file (legacy or no-legacy)")
    parser.add_argument("--snpdists",  required=True, help="Path to snp-dists output file")
    parser.add_argument("--output",    default="mst_visualizer.html", help="Output HTML file (default: mst_visualizer.html)")
    parser.add_argument("--title",     default="Bacterial MST Visualizer", help="Page / report title")
    parser.add_argument("--threshold", type=int, default=None,
                        help="Initial SNP threshold for clustering (default: 5%% of max distance)")
    args = parser.parse_args()

    print(f"[1/4] Parsing MLST file:      {args.mlst}")
    mlst_data = parse_mlst(args.mlst)
    print(f"      → {len(mlst_data)} samples parsed")
    no_st = sum(1 for v in mlst_data.values() if v["st"] == "-")
    if no_st:
        print(f"      ⚠  {no_st} sample(s) have no ST (will be shown in grey)")

    print(f"[2/4] Parsing snp-dists file: {args.snpdists}")
    snp_data = parse_snpdists(args.snpdists)
    print(f"      → {len(snp_data['names'])} samples in matrix")

    max_snp = get_max_snp(snp_data)
    default_thresh = args.threshold if args.threshold is not None else max(1, round(max_snp * 0.05))
    print(f"      Max SNP distance: {max_snp:,}  |  Default threshold: {default_thresh:,}")

    # Warn about samples in SNP matrix but not in MLST
    snp_only = [n for n in snp_data["names"] if n not in mlst_data]
    if snp_only:
        print(f"      ⚠  Samples in SNP matrix but not in MLST: {', '.join(snp_only)}")
        print(f"         These will be shown as grey / no-ST nodes.")

    print(f"[3/4] Building HTML...")
    html = HTML_TEMPLATE.format(
        page_title=args.title,
        mlst_json=json.dumps(mlst_data, ensure_ascii=False),
        snp_names_json=json.dumps(snp_data["names"], ensure_ascii=False),
        snp_matrix_json=json.dumps(snp_data["matrix"], ensure_ascii=False),
        max_snp=max_snp,
        default_thresh=default_thresh,
    )

    out_path = Path(args.output)
    out_path.write_text(html, encoding="utf-8")
    print(f"[4/4] Saved → {out_path.resolve()}")
    print(f"\nDone! Open {args.output} in any modern browser.")


if __name__ == "__main__":
    main()
