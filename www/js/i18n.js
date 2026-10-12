/* Lightweight application-wide language switcher for static and Shiny-rendered text. */
(function () {
  var translations = {
    ru: {
      "Language": "Язык",
      "English": "English",
      "O2T ClinReport - Drug Prioritization Explorer": "O2T ClinReport - исследователь приоритизации лекарств",
      "Drug network - pharmacology - safety - molecular properties - target genes": "лекарственная сеть - фармакология - безопасность - молекулярные свойства - гены-мишени",
      "PDF Report": "PDF-отчёт",
      "Build the drug list and adjust it manually when needed.": "Сформируйте список лекарств и при необходимости скорректируйте его вручную.",
      "Report title": "Название отчёта",
      "Number of drugs (Top N)": "Количество лекарств (Top N)",
      "Network graphs (Top K)": "Сетевых графиков (Top K)",
      "Expert (optional)": "Эксперт (необязательно)",
      "Expert comment (optional)": "Комментарий эксперта (необязательно)",
      "For example: excluded because of a contraindication.": "Например: препарат исключён из-за противопоказаний.",
      "Reset to automatic list": "Сбросить к автоматическому списку",
      "Generate PDF": "Сформировать PDF",
      "Drug list composition and order": "Состав и порядок списка лекарств",
      "Include or exclude drugs and change their order with the arrows. Unchecked drugs remain available as candidates.": "Включайте и исключайте препараты и меняйте их порядок стрелками. Неотмеченные препараты останутся среди кандидатов.",
      "Rank": "Ранг",
      "In PDF": "В PDF",
      "Drug": "Лекарство",
      "Score": "Оценка",
      "Status": "Статус",
      "Order": "Порядок",
      "Expert changes were made to the drug list": "В список лекарств внесены экспертские изменения",
      "No expert changes were made": "Экспертская корректировка не внесена",
      "order changed": "изменён порядок",
      "added:": "добавлены:",
      "removed:": "исключены:",
      "Automatic list for selected Top": "Автоматический список для выбранного Top",
      "At least one drug must remain in the drug list.": "В списке лекарств должно остаться хотя бы одно лекарство.",
      "Analysis inputs": "Входные данные анализа",
      "Patient ID": "ID пациента",
      "Select a drug to inspect the evidence available for its targets, expression profile and interaction network.": "Выберите препарат, чтобы посмотреть данные по его мишеням, профилю экспрессии и сети взаимодействий.",
      "Select a drug to inspect its details": "Выберите препарат, чтобы посмотреть подробности",
      "Interpretation of drug selection": "Интерпретация выбора лекарств",
      "The recommended medicines were prioritised using differential expression": "Рекомендуемые лекарства ранжированы с использованием дифференциальной экспрессии",
      "Pathway enrichment (GSEA) provides metabolic and functional context, while variant-level exome results are reported as supplementary evidence.": "Обогащение биологических путей (GSEA) даёт метаболический и функциональный контекст, а результаты экзомного анализа представлены как дополнительные данные.",
      "Analyses used to build the Drug module": "Анализы, использованные для построения модуля лекарств",
      "Additional analyses and context": "Дополнительные анализы и контекст",
      "Differential expression": "Дифференциальная экспрессия",
      "Used": "Используется",
      "Available": "Доступно",
      "Additional": "Дополнительно",
      "TF regulatory context": "Регуляторный контекст TF",
      "WES variants": "Варианты WES",
      "GSEA analysis": "Анализ GSEA",
      "Not used": "Не используется",
      "Shown as supplementary regulatory context": "Показано как дополнительный регуляторный контекст",
      "Shown as supplementary variant evidence": "Показано как дополнительные данные по вариантам",
      "Network expansion was disabled for this report": "Расширение сети отключено для этого отчёта",
      "Analysis metadata": "Метаданные анализа",
      "Data sources": "Источники данных",
      "Drug annotations": "Аннотации лекарств",
      "Not supplied": "Не предоставлено",
      "Analysis settings": "Параметры анализа",
      "Preparation command": "Команда подготовки",
      "Not recorded": "Не записана",
      "Drug–gene network": "Лекарственно-генная сеть",
      "Drug-gene interaction network": "Сеть взаимодействий лекарство-ген",
      "Drug network": "Лекарственная сеть",
      "records · select a drug to inspect details": "записей · выберите препарат для просмотра подробностей",
      "Drug module": "Модуль лекарств",
      "Network legend": "Легенда сети",
      "Drug information": "Информация о лекарстве",
      "Molecular evidence": "Молекулярные данные",
      "Drug overview": "Обзор лекарства",
      "Mechanism & indications": "Механизм действия и показания",
      "Safety": "Безопасность",
      "Direct targets": "Прямые мишени",
      "WES evidence": "Данные WES",
      "Analysis metadata": "Метаданные анализа",
      "Data sources": "Источники данных",
      "Analysis settings": "Параметры анализа",
      "Preparation command": "Команда подготовки",
      "Show genes": "Показать гены",
      "Expression direction": "Направление экспрессии",
      "Target-linked genes": "Гены, связанные с мишенями",
      "Context genes": "Контекстные гены",
      "Select all": "Выбрать все",
      "Clear": "Очистить",
      "All target-linked pathways": "Все пути, связанные с мишенями",
      "No direct drug targets in this dataset.": "В этом наборе данных нет прямых мишеней лекарств.",
      "Detailed information is not available": "Подробная информация недоступна",
      "for this drug.": "для этого препарата.",
      "Gene profile": "Профиль гена",
      "No genes passed the current filter.": "Ни один ген не прошёл текущий фильтр.",
      "No data available": "Данные недоступны",
      "Collapse drug list": "Свернуть список лекарств"
    }
  };

  var originalText = new WeakMap();

  function translateTextNode(node, language) {
    if (!originalText.has(node)) originalText.set(node, node.nodeValue);
    var source = originalText.get(node);
    if (language === "en") {
      node.nodeValue = source;
      return;
    }
    var value = translations[language] && translations[language][source.trim()];
    if (value) {
      var leading = source.match(/^\s*/)[0];
      var trailing = source.match(/\s*$/)[0];
      node.nodeValue = leading + value + trailing;
      return;
    }
    var dynamic = source.trim();
    if (language === "ru" && dynamic.indexOf("Automatic list for selected Top") === 0) {
      node.nodeValue = source.replace("Automatic list for selected Top", "Автоматический список для выбранного Top").replace("is ready to generate.", "готов к формированию.");
    } else if (language === "ru" && dynamic.indexOf("Expert changes were made to the drug list") === 0) {
      node.nodeValue = source.replace("Expert changes were made to the drug list", "В список лекарств внесены экспертские изменения");
    } else if (language === "ru") {
      var translated = source
        .replace("records · select a drug to inspect details", "записей · выберите препарат для просмотра подробностей")
        .replace(" drug records", " записей лекарств")
        .replace(" genes", " генов")
        .replace("high-confidence edges", "высоконадежных связей")
        .replace("tested gene sets provide pathway context for targets", "протестированных наборов генов дают контекст путей для мишеней")
        .replace("PASS variants across", "PASS вариантов в")
        .replace("WES-specific views are hidden", "специальные представления WES скрыты")
        .replace("Network: ", "Сеть: ")
        .replace(" neighbours: ", " соседей: ")
        .replace("Drug targets, differential expression, and selected pathway coverage", "Мишени лекарства, дифференциальная экспрессия и покрытие выбранных путей")
        .replace("Select a drug in the table to display its", "Выберите препарат в таблице, чтобы показать его")
        .replace("subgraph.", "подграф.");
      if (translated !== source) node.nodeValue = translated;
    }
  }

  function translate(root, language) {
    if (!root) return;
    var walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
    var node;
    while ((node = walker.nextNode())) {
      if (node.parentElement && ["SCRIPT", "STYLE", "OPTION"].indexOf(node.parentElement.tagName) !== -1) continue;
      translateTextNode(node, language);
    }
    if (root.querySelectorAll) {
      root.querySelectorAll("[title]").forEach(function (element) {
        if (!element.dataset.i18nTitle) element.dataset.i18nTitle = element.title;
        if (language === "ru" && translations.ru[element.dataset.i18nTitle]) element.title = translations.ru[element.dataset.i18nTitle];
        if (language === "en") element.title = element.dataset.i18nTitle;
      });
      root.querySelectorAll("[placeholder]").forEach(function (element) {
        if (!element.dataset.i18nPlaceholder) element.dataset.i18nPlaceholder = element.placeholder;
        if (language === "ru" && translations.ru[element.dataset.i18nPlaceholder]) {
          element.placeholder = translations.ru[element.dataset.i18nPlaceholder];
        }
        if (language === "en") element.placeholder = element.dataset.i18nPlaceholder;
      });
    }
  }

  function applyLanguage() {
    var input = document.getElementById("language");
    var selectize = input && input.selectize;
    var language = selectize ? selectize.getValue() : (input ? input.value : "en");
    translate(document.body, language || "en");
    document.documentElement.setAttribute("lang", language || "en");
  }

  document.addEventListener("change", function (event) {
    if (event.target && event.target.id === "language") applyLanguage();
  });
  if (window.jQuery) {
    $(document).on("shiny:connected shiny:value draw.dt", function () {
      window.setTimeout(applyLanguage, 0);
    });
    $(document).on("shiny:inputchanged", function (event, data) {
      if (data && data.name === "language") window.setTimeout(applyLanguage, 0);
    });
  }

  function start() {
    applyLanguage();
    if (window.MutationObserver && document.body) {
      var scheduled = false;
      new MutationObserver(function () {
        if (scheduled) return;
        scheduled = true;
        window.setTimeout(function () {
          scheduled = false;
          applyLanguage();
        }, 0);
      }).observe(document.body, { childList: true, subtree: true });
    }
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
  else start();
})();
