// Open the links of reference lists (DOIs, URLs) in a new tab.
// Reference lists are produced by pandoc from the bibliography of each
// article: entries have the class "csl-entry".
document.addEventListener("DOMContentLoaded", function () {
  document.querySelectorAll(".csl-entry a[href]").forEach(function (link) {
    link.setAttribute("target", "_blank");
    link.setAttribute("rel", "noopener noreferrer");
  });
});
