/* Sugi Variant — set page: client-side sort, filter, flagged-only, copy-link.
   Everything the page needs is already in the DOM; this only reorders/hides rows. */
(function () {
  // ── copy the permalink (the URL *is* the set) ──────────────────────────────
  var copy = document.getElementById("setCopy");
  if (copy) copy.addEventListener("click", function () {
    var restore = copy.textContent;
    function done() {
      copy.textContent = "✓ Copied";
      copy.classList.add("is-copied");
      setTimeout(function () { copy.textContent = restore; copy.classList.remove("is-copied"); }, 1500);
    }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(location.href).then(done, done);
    } else {
      var ta = document.createElement("textarea");
      ta.value = location.href; document.body.appendChild(ta); ta.select();
      try { document.execCommand("copy"); } catch (e) {}
      document.body.removeChild(ta); done();
    }
  });

  var table = document.getElementById("setTable");
  if (!table) return;
  var tbody = table.tBodies[0];
  var rows = Array.prototype.slice.call(tbody.rows);
  var ths = table.tHead.rows[0].cells;

  // ── click a header to sort by that column (toggles asc/desc) ───────────────
  var cur = { col: -1, dir: 1 };
  Array.prototype.forEach.call(ths, function (th, ci) {
    if (!th.dataset.sort) return;
    th.classList.add("sortable");
    th.addEventListener("click", function () {
      var numeric = th.dataset.sort === "num";
      cur = { col: ci, dir: cur.col === ci ? -cur.dir : 1 };
      rows.slice().sort(function (a, b) {
        var av = a.cells[ci].dataset.sv, bv = b.cells[ci].dataset.sv;
        if (av == null) av = a.cells[ci].textContent.trim();
        if (bv == null) bv = b.cells[ci].textContent.trim();
        if (numeric) return ((parseFloat(av) || 0) - (parseFloat(bv) || 0)) * cur.dir;
        return String(av).localeCompare(String(bv)) * cur.dir;
      }).forEach(function (r) { tbody.appendChild(r); });          // reflow in new order
      Array.prototype.forEach.call(ths, function (h) { h.removeAttribute("data-dir"); });
      th.setAttribute("data-dir", cur.dir > 0 ? "asc" : "desc");
    });
  });

  // ── text filter + flagged-only toggle ──────────────────────────────────────
  var filter = document.getElementById("setFilter");
  var flaggedOnly = document.getElementById("setFlaggedOnly");
  var noMatch = document.getElementById("setNoMatch");
  function apply() {
    var q = (filter && filter.value || "").trim().toLowerCase();
    var fo = flaggedOnly && flaggedOnly.checked;
    var shown = 0;
    rows.forEach(function (r) {
      var vis = (!q || r.textContent.toLowerCase().indexOf(q) !== -1) &&
                (!fo || r.dataset.flag === "1");
      r.hidden = !vis;
      if (vis) shown++;
    });
    if (noMatch) noMatch.hidden = shown !== 0;
  }
  if (filter) filter.addEventListener("input", apply);
  if (flaggedOnly) flaggedOnly.addEventListener("change", apply);
})();
