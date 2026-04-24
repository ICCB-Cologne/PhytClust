/* ============================================================
   PhytClust – ui/collapsibles.js
   Collapsible section helper.
   ============================================================ */

export function wireCollapse(toggleId, bodyId) {
  var toggle = document.getElementById(toggleId);
  var body = document.getElementById(bodyId);
  if (!toggle || !body) return;
  toggle.addEventListener("click", function () {
    var isOpen = body.classList.contains("open");
    body.classList.toggle("open", !isOpen);
    toggle.setAttribute("aria-expanded", !isOpen);
    var chev = toggle.querySelector(".chevron");
    if (chev) chev.innerHTML = isOpen ? "&#x25B6;" : "&#x25BC;";
  });
}
