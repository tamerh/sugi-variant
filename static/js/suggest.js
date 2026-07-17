/* Sugi Variant — search autocomplete. Attaches to every search box; queries
   /suggest.json (genes on a single token, variants on GENE + partial change).

   Interaction model:
     • Enter (no suggestion highlighted) → submit the box: resolves 1 item to
       its variant page, or a comma-separated list to a /set page.
     • Pick a VARIANT suggestion (click / Enter on highlight) → adds it to the box
       as a comma-separated token so you can keep building a set, then Enter to go.
     • Pick a GENE suggestion → navigates straight to the gene hub. */
(function () {
  var BASE = window.BASE || "";
  var box = document.createElement("div");
  box.className = "ac-box";
  box.setAttribute("role", "listbox");
  box.style.display = "none";
  document.body.appendChild(box);

  var items = [], sel = -1, curInput = null, timer = null;

  function hide() { box.style.display = "none"; sel = -1; items = []; }
  function esc(s) {
    return (s || "").replace(/[&<>"]/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c];
    });
  }
  // The dropdown aligns to the visible search box, not the (inset) <input> —
  // the hero box pads its input in by ~3.6rem, so anchoring to the input made
  // the menu narrower than and offset from the box.
  function anchorOf(input) {
    return input.closest(".hero-search-box, .header-search") || input;
  }
  function position(input) {
    var r = anchorOf(input).getBoundingClientRect();
    box.style.left = r.left + "px";
    box.style.top = (r.bottom + 5) + "px";
    box.style.width = Math.max(r.width, 240) + "px";
  }
  // The box may hold a comma-separated list; suggestions are for the token
  // currently being typed (everything after the last comma).
  // Grow the hero <textarea> to fit its content (CSS caps it at ~4 lines, then scrolls).
  function autogrow(el) {
    if (el.tagName !== "TEXTAREA") return;
    el.style.height = "auto";
    el.style.height = el.scrollHeight + "px";
  }
  function tokens(input) { return input.value.split(","); }
  function currentToken(input) {
    var t = tokens(input);
    return t[t.length - 1].trim();
  }
  function render() {
    if (!items.length) { hide(); return; }
    box.innerHTML = items.map(function (it, i) {
      var isGene = it.kind === "gene";
      var cls = "ac-item" + (i === sel ? " ac-sel" : "") + (isGene ? "" : " ac-add");
      return '<a class="' + cls + '" role="option" data-i="' + i +
        '" href="' + BASE + "/" + it.url + '">' +
        '<span class="ac-kind ac-' + it.kind + '">' + it.kind + "</span>" +
        '<span class="ac-label">' + esc(it.label) + "</span>" +
        (it.sub ? '<span class="ac-sub">' + esc(it.sub) + "</span>" : "") +
        '<span class="ac-act">' + (isGene ? "Open" : "Add") + "</span></a>";
    }).join("");
    box.style.display = "block";
  }
  function query(input) {
    var q = currentToken(input);
    if (q.length < 2) { hide(); return; }
    fetch(BASE + "/suggest.json?q=" + encodeURIComponent(q))
      .then(function (r) { return r.json(); })
      .then(function (d) {
        if (input !== curInput) return;
        items = d || []; sel = -1; position(input); render();
      }).catch(hide);
  }
  // Pick a suggestion: genes navigate; variants append to the box for a set.
  function pick(input, i) {
    var it = items[i];
    if (!it) return;
    if (it.kind === "gene") { window.location.href = BASE + "/" + it.url; return; }
    var t = tokens(input).map(function (s) { return s.trim(); });
    t.pop();                                   // drop the partial being typed
    t = t.filter(Boolean);
    t.push(it.label);
    input.value = t.join(", ") + ", ";
    autogrow(input);
    hide();
    input.focus();
    position(input);
  }
  function onInput(e) {
    curInput = e.target; autogrow(curInput); position(curInput);
    clearTimeout(timer);
    timer = setTimeout(function () { query(curInput); }, 110);
  }
  function onKey(e) {
    var open = box.style.display !== "none" && items.length;
    if (e.key === "ArrowDown") { if (!open) return; e.preventDefault(); sel = Math.min(sel + 1, items.length - 1); render(); }
    else if (e.key === "ArrowUp") { if (!open) return; e.preventDefault(); sel = Math.max(sel - 1, -1); render(); }
    else if (e.key === "Enter") {
      if (open && sel >= 0) { e.preventDefault(); pick(e.target, sel); }  // add / open the highlighted item
      else if (e.target.tagName === "TEXTAREA") {                         // submit instead of a newline
        e.preventDefault();
        var f = e.target.form;
        if (f) { if (f.requestSubmit) f.requestSubmit(); else f.submit(); }
      }
      // (a plain <input> submits its form on Enter natively)
    }
    else if (e.key === "Escape") { hide(); }
  }
  Array.prototype.forEach.call(document.querySelectorAll('input[type=search][name=q], textarea.hero-search-input'), function (inp) {
    inp.setAttribute("autocomplete", "off");
    inp.addEventListener("input", onInput);
    inp.addEventListener("keydown", onKey);
    inp.addEventListener("blur", function () { setTimeout(hide, 150); });
    inp.addEventListener("focus", function (e) { curInput = e.target; if (currentToken(inp).length >= 2) query(inp); });
  });
  box.addEventListener("mousedown", function (e) { e.preventDefault(); });   // keep input focus on click
  box.addEventListener("click", function (e) {
    var a = e.target.closest(".ac-item");
    if (!a || !curInput) return;
    e.preventDefault();
    pick(curInput, +a.getAttribute("data-i"));
  });
  window.addEventListener("resize", function () { if (curInput) position(curInput); });
  window.addEventListener("scroll", function () { if (curInput && box.style.display !== "none") position(curInput); }, true);
})();
