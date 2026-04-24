# pdb

Ранняя toolkit-заготовка для agentic workflows, которые извлекают и нормализуют данные о `receptor_structure`, `binding_site`, `ligand_complex` и `ligand` из `rcsb.org` и связанных registry.

Проект нужно воспринимать так:

- в первую очередь это data-extraction project
- во вторую очередь это schema-driven normalization project
- в будущем это может стать agent platform, но сейчас это еще не полноценная система

Главная цель репозитория: помогать человеку и coding agent получать строгие, повторяемые, schema-aligned JSON records для основных сущностей проекта.

Главный source of truth для семантики находится в `docs/` и `schemas/`, а не в текущем поведении отдельных script.

## Структура Репозитория

- `README.md` - краткий обзор проекта
- `docs/entity_schema_readme.ru.md` - русское руководство по чтению и применению schema
- `docs/entity_schema_readme.en.md` - английское руководство по чтению и применению schema
- `schemas/*.schema.*.yaml` - schema-файлы сущностей для normalized storage
- `schemas/domain_ontology.master.*.yaml` - ontology index и связи между сущностями
- `tools/` - рабочая зона для extractors, normalization scripts и utility code
- `tools/fetch_receptor_structure.py` - получает одну receptor structure из RCSB и пишет normalized JSON
- `tools/fetch_related_holo_receptor_structures.py` - ищет связанные human/mouse holo structures и переиспользует base fetcher
- `tools/fetch_receptor_by_gene_name.py` - ищет структуры по gene symbol и сохраняет normalized receptor records
- `tools/rcsb_receptor_utils.py` - общие helpers для RCSB Data API, Search API и JSON output
- `output/` - generated JSON outputs на одном уровне с `AGENTS.md`; внутри записи распределяются по папкам `output/<gene_name>/`

`tools/` - это workspace для extraction и normalization logic.

`output/` содержит generated artifacts, а не source code. Каталог находится на одном уровне с `AGENTS.md`, а generated records для отдельных белков должны складываться в gene-scoped subdirectories вроде `output/ACHE/`.

## Сущности Проекта

В проекте сейчас формализованы четыре основные сущности:

- `receptor_structure` - одна normalized запись о receptor/protein structure
- `binding_site` - одна normalized интерпретация residue-level сайта связывания, привязанная к конкретному структурному контексту
- `ligand_complex` - одна normalized интерпретация receptor-ligand complex
- `ligand` - одна normalized chemical ligand entity

Мастер-онтология в `schemas/domain_ontology.master.*.yaml` фиксирует границы этих сущностей и базовые отношения между ними, включая:

- `has_ligand_complex`
- `has_binding_site`
- `has_state_record`
- `derived_from_structure`
- `maps_to_binding_site`
- `has_ligand`
- `participates_in`

## Как Читать Проект

При extraction или normalization работа должна идти в таком порядке:

1. сначала прочитать `docs/entity_schema_readme.ru.md` или `docs/entity_schema_readme.en.md`
2. затем выбрать target schema в `schemas/`
3. воспринимать `fields`, `allowed_values` и `interpretation_rules` как semantic contract
4. формировать JSON, совместимый с `output_template` выбранной schema
5. связывать related records через `related_entities`

Нельзя начинать с догадок о смысле полей только по текущему поведению scripts.

## Правила Нормализации

Ключевые правила, уже зафиксированные в документации и схемах:

- сохранять field names точно как в schema
- сохранять enum-like values точно как в schema
- использовать `null`, если значение отсутствует или не подтверждается trusted source
- не делать silent guessing
- сохранять provenance через registry и identifier fields
- строго сохранять entity boundaries и не смешивать `receptor_structure`, `binding_site`, `ligand_complex` и `ligand`
- связывать сущности через `related_entities`
- относиться консервативно к solvent, salts, buffers и другим non-functional additives

Это особенно важно для интерпретации `ligand_state`, отбора meaningful ligands и различения `apo`/`holo` состояний.

Текущая практическая логика в `tools/` для `apo`/`holo` пока остается несовершенной. Она уже использует nonpolymer entity data из RCSB и умеет находить candidate ligand codes лучше, чем простая проверка `nonpolymer_bound_components`, но все еще может:

- включать в `bound_ligands` сопутствующие non-functional components
- недостаточно надежно выделять один primary ligand code для small-molecule complexes
- ошибаться на границе между `apo`, peptide/protein-bound state и small-molecule holo state

Поэтому текущий extractor нужно считать полезным, но еще не окончательным воплощением schema semantics для `ligand_state`.

Для `receptor_structure` current intent такой:

- `bound_ligands` хранит только meaningful ligand codes после консервативной фильтрации additives
- `bound_ligands` должен по возможности содержать главный ligand среди meaningful candidates, чтобы downstream workflow мог использовать координаты этих лигандов для восстановления protein-ligand binding region после выравнивания структур
- `primary_ligand_code` хранит один функционально центральный ligand code, если его можно выбрать без guessing
- если meaningful ligands несколько и явного winner нет, `primary_ligand_code` должен быть `null`
- при ambiguous cases extractor должен предпочитать `unknown` и explanatory `notes`, а не ложную уверенность

Для `binding_site` current intent такой:

- `residue_members` хранит консервативный набор остатков сайта в конкретном структурном контексте
- `active_residue_members` хранит более строгий поднабор, а не все ligand-contacting residues
- software-generated или просто ligand-contact evidence недостаточно для автоматического вывода о catalytic роли
- если evidence ограничено, интерпретация сайта должна оставаться консервативной, а неопределенность должна явно отражаться в `notes`

## Текущее Состояние

Сейчас это не complete agent pipeline, а ранняя, но уже рабочая основа.

В репозитории уже есть:

- рабочие schema-files
- ontology master
- functional extractor `fetch_receptor_structure.py`
- scripts для поиска связанных holo-структур и структур по gene name
- generated JSON outputs в `output/`

В проекте пока еще нет:

- unified CLI или single entry point для всех workflow
- JSON validation against YAML schema
- extractor или curator workflow для `binding_site`
- extractor для `ligand_complex`
- extractor для `ligand`
- полной normalization coverage в текущих scripts
- dedicated automated test suite

Отдельно: логика разделения `apo`/`holo` и правила выделения key ligand code все еще требуют доработки. Сейчас extractor уже умеет находить более полный набор nonpolymer component codes, но semantic prioritization главного лиганда пока неполная.

Не стоит переоценивать зрелость проекта: текущие scripts полезны, но они еще не покрывают весь schema contract.

## Output Policy

Все generated outputs должны сохраняться в `output/`, который находится на одном уровне с `AGENTS.md`.

Для отдельных белков записи должны раскладываться в gene-scoped поддиректории:

- `output/ACHE/`
- `output/BCHE/`
- `output/CES1/`

Типовой файл receptor structure сейчас выглядит так:

- `output/ACHE/receptor_structure_pdb_5U09.json`

`output/` следует рассматривать как generated artifacts. Эту папку не нужно путать с source code.

## Практические Команды

Так как formal build system пока отсутствует, основной workflow сейчас состоит из запуска scripts и базовой syntax validation.

Проверка окружения:

```bash
python --version
python -m pip --version
```

Получить одну normalized receptor structure:

```bash
python tools/fetch_receptor_structure.py 4EY7 --output-dir output
```

Получить связанные human и mouse holo structures:

```bash
python tools/fetch_related_holo_receptor_structures.py 4EY7 --output-dir output
```

Получить receptor structures по gene symbol:

```bash
python tools/fetch_receptor_by_gene_name.py ACHE --output-dir output
```

Проверить syntax для всего `tools/`:

```bash
python -m compileall tools
```

Проверить отдельные scripts:

```bash
python -m py_compile tools/fetch_receptor_structure.py
python -m py_compile tools/fetch_related_holo_receptor_structures.py
python -m py_compile tools/fetch_receptor_by_gene_name.py
python -m py_compile tools/rcsb_receptor_utils.py
```

## Что Полезно Делать Дальше

Наиболее естественные следующие шаги:

- добавить extractor для `ligand_complex` в `tools/`
- добавить extractor для `ligand` в `tools/`
- добавить extractor или curator workflow для `binding_site`
- добавить JSON validation against schema expectations
- улучшить holo vs apo discrimination
- улучшить выделение primary ligand code для receptor structures
- улучшить normalization coverage в текущем `receptor_structure` extractor
- добавить summary или index outputs для bulk runs
- ввести unified CLI, когда несколько extractors стабилизируются
- добавить lightweight automated tests для известных PDB examples
