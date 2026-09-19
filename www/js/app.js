/* Cross-widget highlights for the selected network gene. */
Shiny.addCustomMessageHandler("filter-drug-rows", function (message) {
  var tableElement = $("#drug_table");
  var allowedRows = message.rows;

  var applyFilter = function () {
    var table = tableElement.DataTable();
    table.settings()[0].geneFilterRows = allowedRows;
    table.draw();
  };

  if ($.fn.dataTable.isDataTable(tableElement)) {
    applyFilter();
  }
});

$.fn.dataTable.ext.search.push(function (settings, data, dataIndex) {
  if (settings.nTable.id !== "drug_table" || settings.geneFilterRows == null) {
    return true;
  }

  return settings.geneFilterRows.indexOf(dataIndex) !== -1;
});

/* Keep the gene-focused drug results adjacent to the DataTable search control. */
(function () {
  function placeFocusedDrugs() {
    var focusedDrugs = document.getElementById("focused_drugs");
    var search = document.getElementById("drug_table_filter");

    if (!focusedDrugs || !search || focusedDrugs.previousElementSibling === search) {
      return;
    }

    search.insertAdjacentElement("afterend", focusedDrugs);
  }

  $(document).on("shiny:connected", placeFocusedDrugs);
  $(document).on("draw.dt", "#drug_table", placeFocusedDrugs);
  $(document).on("shiny:value", "#focused_drugs", placeFocusedDrugs);
})();
