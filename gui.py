#!/usr/bin/env python3
"""
Simplex Method with Bland's Rule - Modern Web GUI
=================================================
A local graphical web interface running on Python standard library.
Opens automatically in your default browser.

Usage:
  python3 gui.py
  python3 gui.py tableau.csv
"""

import os
import sys
import json
import glob
import socket
import urllib.parse
import webbrowser
from http.server import HTTPServer, BaseHTTPRequestHandler
from fractions import Fraction

from simplex import (
    SimplexTableau,
    LinearProgram,
    load_tableau_from_csv,
    record_tableau_steps,
    record_lp_steps,
    parse_fraction,
    format_fraction,
    format_approx,
    get_display_obj_val,
)

HTML_PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Simplex Method with Bland's Rule</title>
  <style>
    :root {
      --bg: #f8fafc;
      --card-bg: #ffffff;
      --border: #e2e8f0;
      --text: #0f172a;
      --text-muted: #64748b;
      --primary: #2563eb;
      --primary-hover: #1d4ed8;
      --primary-light: #eff6ff;
      --success: #16a34a;
      --success-light: #dcfce7;
      --warning: #d97706;
      --warning-light: #fef3c7;
      --purple: #7c3aed;
      --purple-light: #f5f3ff;
      --pivot-cell: #fde047;
      --pivot-col: #dbeafe;
      --pivot-row: #f3e8ff;
      --danger: #dc2626;
      --danger-light: #fee2e2;
    }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background: var(--bg);
      color: var(--text);
      line-height: 1.5;
      padding-bottom: 60px;
    }
    .header {
      background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
      color: white;
      padding: 24px 32px;
      box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1);
    }
    .header h1 { font-size: 24px; font-weight: 700; margin-bottom: 4px; }
    .header p { color: #94a3b8; font-size: 14px; }
    .container { max-width: 1280px; margin: 0 auto; padding: 24px; }

    /* Tabs */
    .tabs {
      display: flex;
      gap: 8px;
      border-bottom: 2px solid var(--border);
      margin-bottom: 24px;
    }
    .tab-btn {
      padding: 10px 18px;
      font-size: 14px;
      font-weight: 600;
      color: var(--text-muted);
      background: transparent;
      border: none;
      border-bottom: 3px solid transparent;
      cursor: pointer;
      display: flex;
      align-items: center;
      gap: 8px;
      transition: all 0.2s;
    }
    .tab-btn:hover { color: var(--primary); }
    .tab-btn.active {
      color: var(--primary);
      border-bottom-color: var(--primary);
    }
    .tab-content { display: none; }
    .tab-content.active { display: block; }

    /* Cards */
    .card {
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 20px;
      margin-bottom: 20px;
      box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    .card-title {
      font-size: 16px;
      font-weight: 700;
      color: var(--text);
      margin-bottom: 14px;
      display: flex;
      align-items: center;
      justify-content: space-between;
    }

    /* Controls Bar */
    .stepper-bar {
      display: flex;
      align-items: center;
      justify-content: space-between;
      flex-wrap: wrap;
      gap: 12px;
      background: #f1f5f9;
      border-radius: 10px;
      padding: 12px 16px;
      margin-bottom: 18px;
    }
    .btn-group { display: flex; align-items: center; gap: 6px; }
    .btn {
      background: var(--primary);
      color: white;
      border: none;
      padding: 8px 14px;
      border-radius: 6px;
      font-size: 13px;
      font-weight: 600;
      cursor: pointer;
      display: inline-flex;
      align-items: center;
      gap: 6px;
      transition: background 0.15s;
    }
    .btn:hover { background: var(--primary-hover); }
    .btn-secondary {
      background: white;
      color: var(--text);
      border: 1px solid var(--border);
    }
    .btn-secondary:hover { background: #f8fafc; }
    .btn:disabled { opacity: 0.5; cursor: not-allowed; }
    .badge {
      padding: 4px 10px;
      border-radius: 20px;
      font-size: 12px;
      font-weight: 700;
    }
    .badge-blue { background: var(--primary-light); color: var(--primary); }
    .badge-green { background: var(--success-light); color: var(--success); }
    .badge-amber { background: var(--warning-light); color: #92400e; border: 1px solid #f59e0b; }
    .badge-purple { background: var(--purple-light); color: #5b21b6; border: 1px solid #8b5cf6; }
    .badge-red { background: var(--danger-light); color: var(--danger); }

    /* Phase Banners */
    .phase-banner {
      border-radius: 8px;
      padding: 12px 18px;
      margin-bottom: 14px;
      display: flex;
      align-items: center;
      gap: 14px;
      font-size: 13px;
    }
    .phase-banner-p1 {
      background: #fffbeb;
      border: 1.5px solid #f59e0b;
      color: #92400e;
    }
    .phase-banner-p2 {
      background: #faf5ff;
      border: 1.5px solid #8b5cf6;
      color: #5b21b6;
    }
    .phase-banner-std {
      background: #eff6ff;
      border: 1.5px solid #3b82f6;
      color: #1e40af;
    }
    .phase-badge-p1 {
      background: #d97706;
      color: white;
      font-weight: 800;
      font-size: 11px;
      letter-spacing: 0.5px;
      padding: 4px 10px;
      border-radius: 4px;
      text-transform: uppercase;
      white-space: nowrap;
    }
    .phase-badge-p2 {
      background: #7c3aed;
      color: white;
      font-weight: 800;
      font-size: 11px;
      letter-spacing: 0.5px;
      padding: 4px 10px;
      border-radius: 4px;
      text-transform: uppercase;
      white-space: nowrap;
    }
    .phase-badge-std {
      background: #2563eb;
      color: white;
      font-weight: 800;
      font-size: 11px;
      letter-spacing: 0.5px;
      padding: 4px 10px;
      border-radius: 4px;
      text-transform: uppercase;
      white-space: nowrap;
    }

    /* Tableau Table */
    .table-container {
      overflow-x: auto;
      margin-bottom: 16px;
      border: 1px solid var(--border);
      border-radius: 8px;
    }
    table.tableau {
      width: 100%;
      border-collapse: collapse;
      font-size: 13px;
      text-align: right;
    }
    table.tableau th, table.tableau td {
      padding: 10px 14px;
      border: 1px solid var(--border);
      white-space: nowrap;
    }
    table.tableau th {
      background: #f8fafc;
      font-weight: 700;
      color: var(--text);
      text-align: center;
    }
    table.tableau th.rhs-col-header, table.tableau td.rhs-cell {
      background: #f1f5f9;
      font-weight: 700;
      border-right: 2px solid #cbd5e1;
      text-align: right;
    }
    table.tableau td.basis-cell {
      text-align: center;
      font-weight: 700;
      background: #f8fafc;
    }

    /* Cost Row (Row 0) Distinct Coloring */
    table.tableau tr.cost-row-p1 {
      background: #fffbeb !important;
      font-weight: 700;
      border-bottom: 2px solid #f59e0b;
    }
    table.tableau tr.cost-row-p1 td {
      color: #92400e;
    }
    table.tableau tr.cost-row-p1 td.basis-cell {
      background: #fef3c7 !important;
      color: #78350f !important;
      font-weight: 800;
    }

    table.tableau tr.cost-row-p2 {
      background: #faf5ff !important;
      font-weight: 700;
      border-bottom: 2px solid #a855f7;
    }
    table.tableau tr.cost-row-p2 td {
      color: #581c87;
    }
    table.tableau tr.cost-row-p2 td.basis-cell {
      background: #f5f3ff !important;
      color: #4c1d95 !important;
      font-weight: 800;
    }

    /* Pivot Highlights */
    table.tableau td.pivot-col-active, table.tableau th.pivot-col-active {
      background: var(--pivot-col) !important;
    }
    table.tableau tr.pivot-row-active td {
      background: var(--pivot-row);
    }
    table.tableau td.pivot-cell-active {
      background: var(--pivot-cell) !important;
      font-weight: 800 !important;
      color: #854d0e !important;
      box-shadow: inset 0 0 0 2px #eab308;
    }

    /* Analysis Box */
    .analysis-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 16px;
      margin-top: 16px;
    }
    .info-card {
      background: #f8fafc;
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 14px;
    }
    .info-card h4 {
      font-size: 13px;
      font-weight: 700;
      margin-bottom: 8px;
      color: var(--text);
      display: flex;
      align-items: center;
      gap: 6px;
    }
    .info-card ul {
      list-style-type: none;
      font-size: 13px;
    }
    .info-card li { margin-bottom: 4px; }

    /* Final Result Banner */
    .result-banner {
      background: var(--success-light);
      border: 1px solid #86efac;
      border-radius: 10px;
      padding: 16px 20px;
      margin-top: 16px;
      display: flex;
      flex-direction: column;
      gap: 8px;
    }
    .result-banner h3 {
      color: #14532d;
      font-size: 16px;
      display: flex;
      align-items: center;
      gap: 8px;
    }

    /* CSV Editor Table */
    .editor-table {
      width: 100%;
      border-collapse: collapse;
      margin-top: 12px;
    }
    .editor-table th, .editor-table td {
      border: 1px solid var(--border);
      padding: 6px 8px;
    }
    .editor-table input {
      width: 100%;
      border: none;
      padding: 6px;
      font-size: 13px;
      font-family: inherit;
      background: transparent;
      outline: none;
    }
    .editor-table input:focus {
      background: #eff6ff;
      border-radius: 4px;
    }
    .select-input {
      padding: 8px 12px;
      border-radius: 6px;
      border: 1px solid var(--border);
      background: white;
      font-size: 14px;
      font-family: inherit;
    }
    .form-group {
      display: flex;
      flex-direction: column;
      gap: 6px;
      margin-bottom: 12px;
    }
    .form-row {
      display: flex;
      gap: 12px;
      align-items: center;
      flex-wrap: wrap;
      margin-bottom: 12px;
    }
  </style>
</head>
<body>

  <header class="header">
    <div class="container" style="padding: 0;">
      <h1>Simplex Method with Bland's Rule</h1>
      <p>Interactive graphical interface with step-by-step tableau visualization and exact fraction arithmetic</p>
    </div>
  </header>

  <div class="container">

    <!-- Tabs Navigation -->
    <div class="tabs">
      <button class="tab-btn active" onclick="switchTab('visualizer')">📊 Step-by-Step Tableau</button>
      <button class="tab-btn" onclick="switchTab('csv')">📁 CSV Files & Tableau Editor</button>
      <button class="tab-btn" onclick="switchTab('manual')">✍️ Manual LP Input</button>
      <button class="tab-btn" onclick="switchTab('demos')">💡 Benchmark & Demo Problems</button>
    </div>

    <!-- TAB 1: VISUALIZER -->
    <div id="tab-visualizer" class="tab-content active">
      <div class="card">
        
        <!-- Controls Bar -->
        <div class="stepper-bar">
          <div class="btn-group">
            <button class="btn btn-secondary" onclick="firstStep()" title="First step">⏮</button>
            <button class="btn btn-secondary" onclick="prevStep()" title="Previous step">◀ Prev</button>
            <span id="step-label" style="font-weight: 700; font-size: 13px; margin: 0 8px;">Tableau 0 / 0</span>
            <button class="btn btn-secondary" onclick="nextStep()" title="Next step">Next ▶</button>
            <button class="btn btn-secondary" onclick="lastStep()" title="Last step">⏭</button>
          </div>

          <div class="btn-group">
            <button id="play-btn" class="btn btn-secondary" onclick="togglePlay()">▶ Auto-play</button>
            <span id="phase-badge" class="badge badge-purple" style="display: none;">Phase</span>
            <span id="status-badge" class="badge badge-blue">In progress</span>
          </div>
        </div>

        <!-- Prominent Phase Banner Container -->
        <div id="phase-banner-container"></div>

        <div style="margin-bottom: 12px; display: flex; justify-content: space-between; align-items: center;">
          <h2 id="tableau-title" style="font-size: 17px; font-weight: 700;">Initial Tableau</h2>
          <div id="tableau-obj" style="font-size: 14px; font-weight: 700; color: var(--purple);">MIN Z = 0</div>
        </div>

        <!-- Main Tableau Table Container -->
        <div class="table-container">
          <table class="tableau" id="tableau-table">
            <!-- Rendered dynamically: Basis | RHS | var_names... (No z column) -->
          </table>
        </div>

        <!-- Current Basic Solution Banner -->
        <div style="background: #f8fafc; border: 1px solid var(--border); border-radius: 8px; padding: 10px 14px; font-size: 13px; margin-bottom: 16px;">
          <strong>Current Basic Solution:</strong>
          <span id="basic-solution-text" style="color: #334155; margin-left: 8px;"></span>
        </div>

        <!-- Bland's Rule Detailed Step Breakdown -->
        <div class="analysis-grid">
          
          <!-- Card 1: Entering Variable -->
          <div class="info-card">
            <h4>🔍 1. Entering Variable Selection (Bland's Rule)</h4>
            <div id="entering-analysis" style="font-size: 13px;"></div>
          </div>

          <!-- Card 2: Ratio Test & Leaving Variable -->
          <div class="info-card">
            <h4>📐 2. Minimum Ratio Test & Leaving Variable (Bland)</h4>
            <div id="leaving-analysis" style="font-size: 13px;"></div>
          </div>

          <!-- Card 3: Pivot Element & Row Operations -->
          <div class="info-card">
            <h4>⚙️ 3. Row Operations & Pivot</h4>
            <div id="pivot-analysis" style="font-size: 13px;"></div>
          </div>

        </div>

        <!-- Final Result Card (shown on optimal) -->
        <div id="final-result-card" class="result-banner" style="display: none;">
          <h3>🎉 Optimal Solution Reached!</h3>
          <div id="final-result-body" style="font-size: 14px;"></div>
        </div>

      </div>
    </div>

    <!-- TAB 2: CSV MANAGER & EDITOR -->
    <div id="tab-csv" class="tab-content">
      <div class="card">
        <div class="card-title">
          <span>Import Tableau from CSV File in Folder</span>
          <button class="btn" onclick="reloadCsvList()">🔄 Refresh List</button>
        </div>

        <div class="form-row">
          <label style="font-weight: 600; font-size: 14px;">Detected CSV Files:</label>
          <select id="csv-select" class="select-input" style="min-width: 240px;" onchange="loadSelectedCsv()">
            <!-- Options populated dynamically -->
          </select>
          <button class="btn" onclick="solveCurrentCsv()">🚀 Solve This CSV</button>
        </div>

        <div style="margin: 16px 0; border-top: 1px solid var(--border); padding-top: 16px;">
          <h4 style="font-size: 14px; font-weight: 700; margin-bottom: 8px;">Or Upload a CSV File from your computer:</h4>
          <input type="file" id="file-upload" accept=".csv" class="select-input" onchange="handleFileUpload(event)">
        </div>

        <!-- Editable CSV Grid -->
        <div style="margin-top: 20px;">
          <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
            <h4 style="font-size: 14px; font-weight: 700;">Tableau Preview & Interactive Editor (Basis, RHS, variables...):</h4>
            <div class="btn-group">
              <button class="btn btn-secondary" onclick="addRow()">+ Add Row</button>
              <button class="btn btn-secondary" onclick="saveCsvToDisk()">💾 Save Changes to CSV</button>
              <button class="btn" onclick="solveEditedTableau()">🚀 Solve Modified Tableau</button>
            </div>
          </div>
          <div class="table-container">
            <table id="editor-table" class="editor-table">
              <!-- Rendered dynamically -->
            </table>
          </div>
        </div>

      </div>
    </div>

    <!-- TAB 3: MANUAL LP INPUT -->
    <div id="tab-manual" class="tab-content">
      <div class="card">
        <div class="card-title">Manual Linear Programming (LP) Input</div>
        <p style="font-size: 13px; color: var(--text-muted); margin-bottom: 16px;">
          Enter the problem specification. The program will automatically construct the initial tableau
          (including Two-Phase Phase 1 with artificial variables y1, y2, ... if &ge; or = constraints are present) and solve it with Bland's Rule.
        </p>

        <div class="form-row">
          <div class="form-group">
            <label style="font-weight: 600; font-size: 13px;">Optimization Goal:</label>
            <select id="manual-obj-type" class="select-input">
              <option value="min" selected>Minimization (MIN Z)</option>
              <option value="max">Maximization (MAX Z)</option>
            </select>
          </div>

          <div class="form-group">
            <label style="font-weight: 600; font-size: 13px;">Decision Variables (n):</label>
            <input type="number" id="manual-n-vars" class="select-input" value="2" min="1" max="10" onchange="buildManualForm()">
          </div>

          <div class="form-group">
            <label style="font-weight: 600; font-size: 13px;">Constraints (m):</label>
            <input type="number" id="manual-n-constraints" class="select-input" value="3" min="1" max="15" onchange="buildManualForm()">
          </div>

          <div class="form-group" style="justify-content: flex-end;">
            <label style="display: flex; align-items: center; gap: 8px; font-weight: 600; font-size: 13px; cursor: pointer; padding-bottom: 8px;">
              <input type="checkbox" id="manual-force-twophase" checked>
              Force Two-Phase Simplex (y artificial basis for all constraints)
            </label>
          </div>
        </div>

        <div id="manual-form-container" style="margin-top: 16px;">
          <!-- Generated inputs -->
        </div>

        <button class="btn" style="margin-top: 20px; padding: 10px 20px; font-size: 14px;" onclick="solveManualLP()">
          🚀 Solve with Simplex Method (Bland's Rule)
        </button>
      </div>
    </div>

    <!-- TAB 4: DEMO PROBLEMS -->
    <div id="tab-demos" class="tab-content">
      <div class="card">
        <div class="card-title">Benchmark & Demonstration Problems (Bland's Rule)</div>
        <div class="analysis-grid">
          
          <div class="info-card">
            <h4>🔄 Beale's Cycling Problem (1955)</h4>
            <p style="font-size: 13px; color: var(--text-muted); margin-bottom: 10px;">
              The quintessential benchmark problem: under Dantzig's largest-coefficient rule, it cycles indefinitely with period 6.
              With <strong>Bland's Rule</strong>, ties are resolved and it converges to the optimum in 7 iterations!
            </p>
            <button class="btn" onclick="loadDemo('beale')">Load & Solve Beale</button>
          </div>

          <div class="info-card">
            <h4>📈 Standard 2D Problem</h4>
            <p style="font-size: 13px; color: var(--text-muted); margin-bottom: 10px;">
              Max 3x1 + 5x2 subject to standard inequality constraints &le;. Clean pedagogical step-by-step example.
            </p>
            <button class="btn" onclick="loadDemo('standard2d')">Load & Solve 2D</button>
          </div>

          <div class="info-card">
            <h4>⚖️ Two-Phase Problem (Constraints =, &ge;, &le;)</h4>
            <p style="font-size: 13px; color: var(--text-muted); margin-bottom: 10px;">
              Min 4x1 + x2 with mixed constraints forcing artificial variables <strong>y1, y2, y3</strong> as initial basis.
              Demonstrates Phase 1 (Min W &rarr; 0) and Phase 2 (Min Z &rarr; 17/5 = 3.4).
            </p>
            <button class="btn" onclick="loadDemo('twophase')">Load & Solve Two-Phase</button>
          </div>

        </div>
      </div>
    </div>

  </div>

  <script>
    let stepsData = [];
    let currentStepIdx = 0;
    let playInterval = null;
    let editorHeaders = [];
    let editorRows = [];

    function switchTab(tabId) {
      document.querySelectorAll('.tab-btn').forEach(btn => btn.classList.remove('active'));
      document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
      document.querySelector(`[onclick="switchTab('${tabId}')"]`).classList.add('active');
      document.getElementById(`tab-${tabId}`).classList.add('active');
    }

    // Load CSV list on init
    async function reloadCsvList() {
      try {
        const res = await fetch('/api/csv-files');
        const files = await res.json();
        const select = document.getElementById('csv-select');
        select.innerHTML = '';
        files.forEach(f => {
          const opt = document.createElement('option');
          opt.value = f;
          opt.textContent = f;
          select.appendChild(opt);
        });
        if (files.includes('tableau.csv')) {
          select.value = 'tableau.csv';
        } else if (files.length > 0) {
          select.value = files[0];
        }
        if (select.value) {
          await loadSelectedCsv();
        }
      } catch (err) {
        console.error('Error fetching CSV list:', err);
      }
    }

    async function loadSelectedCsv() {
      const select = document.getElementById('csv-select');
      if (!select.value) return;
      try {
        const res = await fetch('/api/load-csv?file=' + encodeURIComponent(select.value));
        const data = await res.json();
        if (data.error) {
          alert('Error loading CSV: ' + data.error);
          return;
        }
        editorHeaders = data.headers;
        editorRows = data.rows;
        renderEditorTable();
      } catch (err) {
        alert('Error loading CSV: ' + err);
      }
    }

    function renderEditorTable() {
      const table = document.getElementById('editor-table');
      table.innerHTML = '';
      
      // Header row
      const thead = document.createElement('thead');
      const trH = document.createElement('tr');
      editorHeaders.forEach((h, cIdx) => {
        const th = document.createElement('th');
        th.innerHTML = `<input type="text" value="${h}" onchange="editorHeaders[${cIdx}] = this.value">`;
        trH.appendChild(th);
      });
      thead.appendChild(trH);
      table.appendChild(thead);

      // Body rows
      const tbody = document.createElement('tbody');
      editorRows.forEach((row, rIdx) => {
        const tr = document.createElement('tr');
        row.forEach((cell, cIdx) => {
          const td = document.createElement('td');
          td.innerHTML = `<input type="text" value="${cell}" onchange="editorRows[${rIdx}][${cIdx}] = this.value">`;
          tr.appendChild(td);
        });
        tbody.appendChild(tr);
      });
      table.appendChild(tbody);
    }

    function addRow() {
      const newRow = new Array(editorHeaders.length).fill('0');
      newRow[0] = 's' + editorRows.length;
      newRow[1] = '0'; // RHS column
      editorRows.push(newRow);
      renderEditorTable();
    }

    async function saveCsvToDisk() {
      const filename = document.getElementById('csv-select').value || 'tableau.csv';
      const content = [editorHeaders.join(',')].concat(editorRows.map(r => r.join(','))).join('\\n');
      try {
        await fetch('/api/save-csv', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ filename, content })
        });
        alert('File ' + filename + ' saved successfully!');
      } catch (e) {
        alert('Error saving file: ' + e);
      }
    }

    async function solveCurrentCsv() {
      const select = document.getElementById('csv-select');
      if (!select.value) return;
      try {
        const res = await fetch('/api/solve-file?file=' + encodeURIComponent(select.value));
        const text = await res.text();
        let data;
        try {
          data = JSON.parse(text);
        } catch (e) {
          console.error("Non-JSON response:", text);
          alert("Server response error: " + text.slice(0, 150));
          return;
        }
        if (data.error) {
          alert("Error solving: " + data.error);
          return;
        }
        setSolutionData(data);
        switchTab('visualizer');
      } catch (err) {
        alert('Error solving: ' + err);
      }
    }

    async function solveEditedTableau() {
      try {
        const res = await fetch('/api/solve-tableau', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ headers: editorHeaders, rows: editorRows })
        });
        const text = await res.text();
        let data;
        try {
          data = JSON.parse(text);
        } catch (e) {
          console.error("Non-JSON response:", text);
          alert("Server response error: " + text.slice(0, 150));
          return;
        }
        if (data.error) {
          alert("Error solving: " + data.error);
          return;
        }
        setSolutionData(data);
        switchTab('visualizer');
      } catch (err) {
        alert('Error solving: ' + err);
      }
    }

    function setSolutionData(data) {
      stepsData = data.steps;
      currentStepIdx = 0;
      renderCurrentStep();
    }

    function renderCurrentStep() {
      if (!stepsData || stepsData.length === 0) return;
      const step = stepsData[currentStepIdx];
      const total = stepsData.length;

      document.getElementById('step-label').textContent = `Step ${currentStepIdx + 1} of ${total}`;
      document.getElementById('tableau-title').textContent = step.title;

      const objLabel = (step.obj_type ? step.obj_type.toUpperCase() + ' ' : '') + step.obj_name.toUpperCase();
      document.getElementById('tableau-obj').textContent = `${objLabel} = ${step.obj_val}`;

      // Stepper Bar Badges
      const phaseBadge = document.getElementById('phase-badge');
      if (step.is_phase1) {
        phaseBadge.style.display = 'inline-block';
        phaseBadge.textContent = '⚡ Phase 1 (Min W)';
        phaseBadge.className = 'badge badge-amber';
      } else if (step.is_phase2) {
        phaseBadge.style.display = 'inline-block';
        phaseBadge.textContent = `🎯 Phase 2 (${objLabel})`;
        phaseBadge.className = 'badge badge-purple';
      } else {
        phaseBadge.style.display = 'inline-block';
        phaseBadge.textContent = `📊 Standard Simplex (${objLabel})`;
        phaseBadge.className = 'badge badge-blue';
      }

      // Status badge
      const badge = document.getElementById('status-badge');
      if (step.status === 'OPTIMAL') {
        badge.textContent = 'Optimal Reached';
        badge.className = 'badge badge-green';
      } else if (step.status === 'UNBOUNDED') {
        badge.textContent = 'Unbounded Problem';
        badge.className = 'badge badge-red';
      } else {
        badge.textContent = `Iteration ${step.iteration}`;
        badge.className = 'badge badge-blue';
      }

      // Phase Banner rendering
      const bannerContainer = document.getElementById('phase-banner-container');
      if (step.is_phase1) {
        bannerContainer.innerHTML = `
          <div class="phase-banner phase-banner-p1">
            <span class="phase-badge-p1">⚡ PHASE 1 ACTIVE</span>
            <div>
              <strong>Phase 1: Finding an Initial Feasible Basis</strong> &mdash;
              Minimizing the sum of artificial variables (<code>y1, y2, ...</code>).
              Current goal: Drive <code>W &rarr; 0</code> so all artificial variables leave the basis.
            </div>
          </div>`;
      } else if (step.is_phase2) {
        bannerContainer.innerHTML = `
          <div class="phase-banner phase-banner-p2">
            <span class="phase-badge-p2">🎯 PHASE 2 ACTIVE</span>
            <div>
              <strong>Phase 2: Optimizing Original Objective</strong> &mdash;
              Feasible basis achieved (<code>W = 0</code>)! All artificial variables have been eliminated.
              Now optimizing original objective <code>${objLabel}</code>.
            </div>
          </div>`;
      } else {
        bannerContainer.innerHTML = `
          <div class="phase-banner phase-banner-std">
            <span class="phase-badge-std">STANDARD SIMPLEX</span>
            <div>
              Optimizing original objective <code>${objLabel}</code> using standard Simplex with Bland's Rule.
            </div>
          </div>`;
      }

      // Basic solution text
      const solParts = [];
      for (const [k, v] of Object.entries(step.basic_solution)) {
        solParts.push(`<strong>${k}</strong> = ${v}`);
      }
      document.getElementById('basic-solution-text').innerHTML = solParts.join(', ');

      // Build tableau HTML
      const table = document.getElementById('tableau-table');
      table.innerHTML = '';

      // Determine pivot row and col indices
      const pivotRow = step.pivot ? step.pivot.row : -1;
      const pivotCol = step.pivot ? step.pivot.col : -1;

      // Table Header: Basis | RHS | var_names... (NO z column)
      const thead = document.createElement('thead');
      const trH = document.createElement('tr');
      trH.innerHTML = `<th>Basis</th><th class="rhs-col-header">RHS</th>`;
      step.var_names.forEach((name, cIdx) => {
        const isPivotCol = cIdx === pivotCol;
        const cls = isPivotCol ? 'pivot-col-active' : '';
        trH.innerHTML += `<th class="${cls}">${name}</th>`;
      });
      thead.appendChild(trH);
      table.appendChild(thead);

      // Table Body
      const tbody = document.createElement('tbody');
      step.rows.forEach((row, rIdx) => {
        const tr = document.createElement('tr');
        if (rIdx === 0) {
          tr.className = step.is_phase1 ? 'cost-row-p1' : 'cost-row-p2';
        }
        if (rIdx === pivotRow) tr.classList.add('pivot-row-active');

        // Basis cell and RHS cell (RHS placed immediately next to Basis)
        const basisBadge = rIdx === 0
          ? `<span class="badge ${step.is_phase1 ? 'badge-amber' : 'badge-purple'}">${row.basic}</span>`
          : `<strong>${row.basic}</strong>`;
        tr.innerHTML = `<td class="basis-cell">${basisBadge}</td><td class="rhs-cell">${row.rhs}</td>`;

        // Variable coefficient cells
        row.cells.forEach((cellVal, cIdx) => {
          let cellCls = '';
          if (rIdx === pivotRow && cIdx === pivotCol) {
            cellCls = 'pivot-cell-active';
          } else if (cIdx === pivotCol) {
            cellCls = 'pivot-col-active';
          }
          tr.innerHTML += `<td class="${cellCls}">${cellVal}</td>`;
        });

        tbody.appendChild(tr);
      });
      table.appendChild(tbody);

      // Analysis Cards
      renderAnalysis(step, objLabel);
    }

    function renderAnalysis(step, objLabel) {
      const enteringDiv = document.getElementById('entering-analysis');
      const leavingDiv = document.getElementById('leaving-analysis');
      const pivotDiv = document.getElementById('pivot-analysis');
      const finalCard = document.getElementById('final-result-card');

      if (step.status === 'OPTIMAL') {
        enteringDiv.innerHTML = '<span class="badge badge-green">No negative reduced cost</span><p style="margin-top: 8px;">All coefficients in Row 0 are &ge; 0. The tableau is optimal!</p>';
        leavingDiv.innerHTML = 'No variable leaves the basis.';
        pivotDiv.innerHTML = 'No row operations needed.';
        
        // Show final result card
        finalCard.style.display = 'block';
        let solHtml = `<p><strong>Optimal Value:</strong> ${objLabel}* = <strong>${step.obj_val}</strong></p><ul style="margin-top: 6px; list-style-type: none;">`;
        for (const [v, val] of Object.entries(step.basic_solution)) {
          solHtml += `<li>• ${v} = <strong>${val}</strong> (Basic)</li>`;
        }
        solHtml += '</ul>';
        document.getElementById('final-result-body').innerHTML = solHtml;
        return;
      }

      finalCard.style.display = 'none';

      // 1. Entering Variable
      if (step.candidates && step.candidates.length > 0) {
        let entHtml = '<ul style="margin-bottom: 8px;">';
        step.candidates.forEach(c => {
          const isChosen = step.entering_var && c.name === step.entering_var.name;
          entHtml += `<li>• <strong>${c.name}</strong> (Index ${c.index}) | Coeff: ${c.coeff} ${isChosen ? '<span class="badge badge-blue">Chosen by Bland</span>' : ''}</li>`;
        });
        entHtml += `</ul><p style="font-size: 12px; color: var(--text-muted);">Bland's Rule selects <strong>${step.entering_var.name}</strong> because it has the smallest index among all candidates.</p>`;
        enteringDiv.innerHTML = entHtml;
      } else {
        enteringDiv.innerHTML = 'No candidates.';
      }

      // 2. Leaving Variable & Ratios
      if (step.ratios && step.ratios.length > 0) {
        let rHtml = '<ul style="margin-bottom: 8px;">';
        step.ratios.forEach(r => {
          const isLeaving = step.leaving_var && r.row === step.leaving_var.row;
          if (r.valid) {
            rHtml += `<li>• Row ${r.row} (${r.basic}, Idx ${r.index}): ${r.rhs} / ${r.coeff} = <strong>${r.ratio}</strong> ${isLeaving ? '<span class="badge badge-purple">Leaves Basis</span>' : ''}</li>`;
          } else {
            rHtml += `<li style="color: var(--text-muted);">• Row ${r.row} (${r.basic}): No ratio (${r.reason})</li>`;
          }
        });
        rHtml += '</ul>';
        if (step.bland_tie_note) {
          rHtml += `<div style="background: var(--warning-light); color: #92400e; padding: 6px 10px; border-radius: 6px; font-size: 12px; margin-top: 6px;"><strong>Degeneracy / Tie:</strong> ${step.bland_tie_note}</div>`;
        }
        leavingDiv.innerHTML = rHtml;
      } else {
        leavingDiv.innerHTML = 'No valid ratios.';
      }

      // 3. Pivot & Operations
      if (step.pivot) {
        let pvtHtml = `<p><strong>Pivot Element:</strong> Row ${step.pivot.row}, Column '${step.entering_var.name}' (Value: <strong>${step.pivot.val}</strong>)</p>`;
        pvtHtml += '<ul style="margin-top: 8px; font-family: monospace; font-size: 12px;">';
        step.row_operations.forEach(op => {
          pvtHtml += `<li>• ${op}</li>`;
        });
        pvtHtml += '</ul>';
        pivotDiv.innerHTML = pvtHtml;
      } else {
        pivotDiv.innerHTML = '-';
      }
    }

    function firstStep() { currentStepIdx = 0; renderCurrentStep(); }
    function prevStep() { if (currentStepIdx > 0) { currentStepIdx--; renderCurrentStep(); } }
    function nextStep() { if (currentStepIdx < stepsData.length - 1) { currentStepIdx++; renderCurrentStep(); } }
    function lastStep() { currentStepIdx = stepsData.length - 1; renderCurrentStep(); }

    function togglePlay() {
      const btn = document.getElementById('play-btn');
      if (playInterval) {
        clearInterval(playInterval);
        playInterval = null;
        btn.textContent = '▶ Auto-play';
      } else {
        btn.textContent = '⏸ Pause';
        playInterval = setInterval(() => {
          if (currentStepIdx < stepsData.length - 1) {
            nextStep();
          } else {
            togglePlay();
          }
        }, 1200);
      }
    }

    // Manual Form Builder
    function buildManualForm() {
      const nVars = parseInt(document.getElementById('manual-n-vars').value) || 2;
      const nConst = parseInt(document.getElementById('manual-n-constraints').value) || 3;
      const container = document.getElementById('manual-form-container');
      container.innerHTML = '';

      // Objective Row
      let objHtml = '<div style="margin-bottom: 16px;"><strong>Objective Function Coefficients (Z):</strong><div style="display: flex; gap: 8px; margin-top: 6px; flex-wrap: wrap;">';
      for (let j = 1; j <= nVars; j++) {
        objHtml += `<div style="display: flex; align-items: center; gap: 4px;"><input type="text" id="obj-c-${j}" class="select-input" style="width: 70px;" value="${j === 1 ? '4' : '1'}"> x${j}</div>`;
      }
      objHtml += '</div></div>';

      // Constraints
      let cHtml = '<strong>Constraints:</strong><div style="display: flex; flex-direction: column; gap: 8px; margin-top: 6px;">';
      const sampleDefs = [
        { c: ['3', '1'], rel: '=', rhs: '3' },
        { c: ['4', '3'], rel: '>=', rhs: '6' },
        { c: ['1', '2'], rel: '<=', rhs: '4' }
      ];
      for (let i = 1; i <= nConst; i++) {
        cHtml += `<div style="display: flex; align-items: center; gap: 8px; flex-wrap: wrap;"><span>(${i})</span>`;
        for (let j = 1; j <= nVars; j++) {
          const defaultC = (i <= sampleDefs.length && j <= sampleDefs[i-1].c.length) ? sampleDefs[i-1].c[j-1] : (j === 1 ? i : '1');
          cHtml += `<input type="text" id="c-${i}-${j}" class="select-input" style="width: 70px;" value="${defaultC}"> x${j} ${j < nVars ? '+' : ''}`;
        }
        const defaultRel = (i <= sampleDefs.length) ? sampleDefs[i-1].rel : '<=';
        const defaultRhs = (i <= sampleDefs.length) ? sampleDefs[i-1].rhs : (i * 10);
        cHtml += `
          <select id="rel-${i}" class="select-input" style="width: 70px;">
            <option value="<=" ${defaultRel === '<=' ? 'selected' : ''}>&le;</option>
            <option value=">=" ${defaultRel === '>=' ? 'selected' : ''}>&ge;</option>
            <option value="=" ${defaultRel === '=' ? 'selected' : ''}>=</option>
          </select>
          <input type="text" id="rhs-${i}" class="select-input" style="width: 80px;" value="${defaultRhs}">
        </div>`;
      }
      cHtml += '</div>';

      container.innerHTML = objHtml + cHtml;
    }

    async function solveManualLP() {
      const objType = document.getElementById('manual-obj-type').value;
      const nVars = parseInt(document.getElementById('manual-n-vars').value);
      const nConst = parseInt(document.getElementById('manual-n-constraints').value);
      const forceTwoPhase = document.getElementById('manual-force-twophase').checked;

      const objCoeffs = [];
      for (let j = 1; j <= nVars; j++) {
        objCoeffs.push(document.getElementById(`obj-c-${j}`).value || '0');
      }

      const constraints = [];
      for (let i = 1; i <= nConst; i++) {
        const coeffs = [];
        for (let j = 1; j <= nVars; j++) {
          coeffs.push(document.getElementById(`c-${i}-${j}`).value || '0');
        }
        const rel = document.getElementById(`rel-${i}`).value;
        const rhs = document.getElementById(`rhs-${i}`).value || '0';
        constraints.push({ coeffs, rel, rhs });
      }

      try {
        const res = await fetch('/api/solve-lp', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ obj_type: objType, obj_coeffs: objCoeffs, constraints, force_two_phase: forceTwoPhase })
        });
        const data = await res.json();
        setSolutionData(data);
        switchTab('visualizer');
      } catch (err) {
        alert('Error solving manual LP: ' + err);
      }
    }

    async function loadDemo(demoName) {
      try {
        const res = await fetch('/api/demo?name=' + demoName);
        const data = await res.json();
        setSolutionData(data);
        switchTab('visualizer');
      } catch (err) {
        alert('Error loading demo: ' + err);
      }
    }

    function handleFileUpload(e) {
      const file = e.target.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = async function(evt) {
        const text = evt.target.result;
        try {
          const res = await fetch('/api/upload-csv', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ filename: file.name, content: text })
          });
          const data = await res.json();
          await reloadCsvList();
          document.getElementById('csv-select').value = file.name;
          loadSelectedCsv();
          alert('File ' + file.name + ' uploaded successfully!');
        } catch (err) {
          alert('Error uploading file: ' + err);
        }
      };
      reader.readAsText(file);
    }

    // Auto load on page ready
    window.onload = async function() {
      buildManualForm();
      await reloadCsvList();
      await solveCurrentCsv();
    };
  </script>
</body>
</html>
"""


class SimplexGUIRequestHandler(BaseHTTPRequestHandler):
    def log_message(self, format, *args):
        # Keep terminal log clean
        pass

    def do_GET(self):
        try:
            parsed = urllib.parse.urlparse(self.path)
            path = parsed.path
            query = urllib.parse.parse_qs(parsed.query)

            if path in ("/", "/index.html"):
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.end_headers()
                self.wfile.write(HTML_PAGE.encode("utf-8"))
                return

            elif path == "/api/csv-files":
                csv_files = sorted(glob.glob("*.csv"))
                self.send_json(csv_files)
                return

            elif path == "/api/load-csv":
                filename = query.get("file", ["tableau.csv"])[0]
                if not os.path.exists(filename):
                    self.send_json({"error": f"File '{filename}' not found."}, status_code=404)
                    return
                tab, obj_type = load_tableau_from_csv(filename)
                headers = ["Basis", "RHS"] + tab.var_names
                rows = []
                display_cost = get_display_obj_val(tab.tableau[0][-1], obj_type, tab.obj_name)
                row0 = [tab.obj_name, format_fraction(display_cost)] + [format_fraction(val) for val in tab.tableau[0][:-1]]
                rows.append(row0)
                for r in range(1, tab.num_rows()):
                    b_var = tab.basic_vars[r - 1]
                    row_r = [b_var, format_fraction(tab.tableau[r][-1])] + [format_fraction(val) for val in tab.tableau[r][:-1]]
                    rows.append(row_r)
                self.send_json({"headers": headers, "rows": rows, "obj_type": obj_type})
                return

            elif path == "/api/solve-file":
                filename = query.get("file", ["tableau.csv"])[0]
                if not os.path.exists(filename):
                    self.send_json({"error": f"File '{filename}' not found."}, status_code=404)
                    return
                tab, obj_type = load_tableau_from_csv(filename)
                result = record_tableau_steps(tab, obj_type)
                self.send_json(result)
                return

            elif path == "/api/demo":
                name = query.get("name", ["beale"])[0]
                if name == "beale":
                    from simplex import get_beale_cycling_example
                    lp = get_beale_cycling_example()
                    result = record_lp_steps(lp)
                elif name == "standard2d":
                    from simplex import get_standard_2d_example
                    lp = get_standard_2d_example()
                    result = record_lp_steps(lp)
                elif name == "twophase":
                    from simplex import get_twophase_example
                    lp = get_twophase_example()
                    result = record_lp_steps(lp)
                else:
                    self.send_json({"error": f"Unknown demo '{name}'."}, status_code=400)
                    return
                self.send_json(result)
                return

            self.send_json({"error": "Endpoint not found."}, status_code=404)
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.send_json({"error": str(e)}, status_code=500)

    def do_POST(self):
        try:
            parsed = urllib.parse.urlparse(self.path)
            path = parsed.path
            length = int(self.headers.get("Content-Length", 0))
            body = self.rfile.read(length).decode("utf-8")
            data = json.loads(body) if body else {}

            if path == "/api/save-csv":
                filename = data.get("filename", "tableau.csv")
                content = data.get("content", "")
                with open(filename, "w", encoding="utf-8") as f:
                    f.write(content)
                self.send_json({"success": True})
                return

            elif path == "/api/upload-csv":
                filename = data.get("filename", "uploaded.csv")
                content = data.get("content", "")
                with open(filename, "w", encoding="utf-8") as f:
                    f.write(content)
                self.send_json({"success": True})
                return

            elif path == "/api/solve-tableau":
                headers = data.get("headers", [])
                rows = data.get("rows", [])
                # If headers[1] is RHS (Basis, RHS, x1, x2, ...)
                if len(headers) > 2 and headers[1].lower() in ('rhs', 'b', 'val', 'costo', 'term'):
                    var_names = headers[2:]
                    matrix = []
                    for r in rows:
                        rhs_val = parse_fraction(r[1])
                        var_vals = [parse_fraction(cell) for cell in r[2:]]
                        matrix.append(var_vals + [rhs_val])
                else:
                    var_names = headers[1:-1]
                    matrix = []
                    for r in rows:
                        matrix.append([parse_fraction(cell) for cell in r[1:]])
                var_indices = {name: idx + 1 for idx, name in enumerate(var_names)}
                obj_name = rows[0][0].strip().lower() or "z"
                basic_vars = [r[0].strip() for r in rows[1:]]
                tab = SimplexTableau(var_names, var_indices, matrix, basic_vars, obj_name=obj_name)
                result = record_tableau_steps(tab, obj_type="min")
                self.send_json(result)
                return

            elif path == "/api/solve-lp":
                obj_type = data.get("obj_type", "min")
                obj_coeffs = [parse_fraction(c) for c in data.get("obj_coeffs", [])]
                raw_constraints = data.get("constraints", [])
                force_two_phase = data.get("force_two_phase", False)
                constraints = []
                for c in raw_constraints:
                    coeffs = [parse_fraction(x) for x in c["coeffs"]]
                    rel = c["rel"]
                    rhs = parse_fraction(c["rhs"])
                    constraints.append((coeffs, rel, rhs))
                lp = LinearProgram(obj_type, obj_coeffs, constraints, force_two_phase=force_two_phase)
                result = record_lp_steps(lp)
                self.send_json(result)
                return

            self.send_json({"error": "Endpoint not found."}, status_code=404)
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.send_json({"error": str(e)}, status_code=500)

    def send_json(self, data, status_code=200):
        self.send_response(status_code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.end_headers()
        def fallback(obj):
            if isinstance(obj, Fraction):
                return format_fraction(obj)
            return str(obj)
        self.wfile.write(json.dumps(data, default=fallback).encode("utf-8"))


def find_free_port(start_port: int = 5050) -> int:
    port = start_port
    while port < 65535:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', port)) != 0:
                return port
        port += 1
    return start_port


def run_gui(port: int = 5050):
    port = find_free_port(port)
    server_address = ('127.0.0.1', port)
    httpd = HTTPServer(server_address, SimplexGUIRequestHandler)
    url = f"http://127.0.0.1:{port}"

    print("=" * 82)
    print("      GRAPHICAL USER INTERFACE (GUI) - SIMPLEX METHOD WITH BLAND'S RULE     ")
    print("=" * 82)
    print(f"Server successfully started at: {url}")
    print("Opening browser automatically...")
    print("Press Ctrl+C in terminal to stop server.")
    print("=" * 82 + "\n")

    try:
        webbrowser.open(url)
    except Exception as e:
        print(f"Could not open browser automatically: {e}")
        print(f"Please open this URL manually: {url}")

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down server...")
        httpd.server_close()
        print("Server stopped.")


if __name__ == "__main__":
    run_gui()
