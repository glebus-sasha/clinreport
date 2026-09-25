/* Cross-widget highlights for the selected network gene. */
$(document).on("click", "#drug_network_toggle", function () {
  var panel = $(this).closest(".embedded-gsea-panel");
  panel.toggleClass("is-collapsed");
  $(this).text(panel.hasClass("is-collapsed") ? "‹" : "›");
});

$(document).on("click", "#drug_list_toggle", function () {
  var layout = $(this).closest(".main-layout");
  layout.toggleClass("drug-list-collapsed");
  $(this).text(layout.hasClass("drug-list-collapsed") ? "›" : "‹");
});

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

/* Select a WES-listed gene in the current interaction graph when it is present. */
Shiny.addCustomMessageHandler("focus-network-gene", function (message) {
  var widget = HTMLWidgets.find("#drug_network_graph");
  if (!widget || !widget.getInstance) return;

  var instance = widget.getInstance();
  var network = instance && instance.network;
  if (!network || !network.body || !network.body.nodes) return;

  var nodeId = network.body.nodes[message.id]
    ? message.id
    : "__PATHWAY_GENE__" + message.symbol;
  if (!network.body.nodes[nodeId]) return;

  network.selectNodes([nodeId], false);
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
