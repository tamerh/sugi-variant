/* Sugi Variant — search autocomplete. Attaches to every search box; queries
   /suggest.json (genes on a single token, variants on GENE + partial change). */
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
  function position(input) {
    var r = input.getBoundingClientRect();
    box.style.left = r.left + "px";
    box.style.top = (r.bottom + 5) + "px";
    box.style.width = Math.max(r.width, 260) + "px";
  }
  function render() {
    if (!items.length) { hide(); return; }
    box.innerHTML = items.map(function (it, i) {
      return '<a class="ac-item' + (i === sel ? " ac-sel" : "") + '" role="option" data-i="' + i +
        '" href="' + BASE + "/" + it.url + '">' +
        '<span class="ac-kind ac-' + it.kind + '">' + it.kind + "</span>" +
        '<span class="ac-label">' + esc(it.label) + "</span>" +
        (it.sub ? '<span class="ac-sub">' + esc(it.sub) + "</span>" : "") + "</a>";
    }).join("");
    box.style.display = "block";
  }
  function query(input) {
    var q = input.value.trim();
    if (q.length < 2) { hide(); return; }
    fetch(BASE + "/suggest.json?q=" + encodeURIComponent(q))
      .then(function (r) { return r.json(); })
      .then(function (d) {
        if (input !== curInput) return;
        items = d || []; sel = -1; position(input); render();
      }).catch(hide);
  }
  function onInput(e) {
    curInput = e.target; position(curInput);
    clearTimeout(timer);
    timer = setTimeout(function () { query(curInput); }, 110);
  }
  function onKey(e) {
    if (box.style.display === "none" || !items.length) return;
    if (e.key === "ArrowDown") { e.preventDefault(); sel = Math.min(sel + 1, items.length - 1); render(); }
    else if (e.key === "ArrowUp") { e.preventDefault(); sel = Math.max(sel - 1, -1); render(); }
    else if (e.key === "Enter" && sel >= 0) { e.preventDefault(); window.location.href = BASE + "/" + items[sel].url; }
    else if (e.key === "Escape") { hide(); }
  }
  Array.prototype.forEach.call(document.querySelectorAll('input[type=search][name=q]'), function (inp) {
    inp.setAttribute("autocomplete", "off");
    inp.addEventListener("input", onInput);
    inp.addEventListener("keydown", onKey);
    inp.addEventListener("blur", function () { setTimeout(hide, 150); });
    inp.addEventListener("focus", function (e) { curInput = e.target; if (inp.value.trim().length >= 2) query(inp); });
  });
  box.addEventListener("mousedown", function (e) { e.preventDefault(); }); // keep input focus on click
  window.addEventListener("resize", function () { if (curInput) position(curInput); });
  window.addEventListener("scroll", function () { if (curInput && box.style.display !== "none") position(curInput); }, true);
})();
