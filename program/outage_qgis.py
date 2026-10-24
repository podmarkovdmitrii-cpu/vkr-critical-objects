# -*- coding: utf-8 -*-
"""
ПРОГРАММА АНАЛИЗА АВАРИЙНЫХ ОТКЛЮЧЕНИЙ ЭЛЕКТРИЧЕСКОЙ СЕТИ  (QGIS 3.x, консоль Python)   версия 2.0
------------------------------------------------------------------------------------
Запуск: Модули -> Консоль Python -> в консоли выполнить:
    exec(open(r'C:\\Users\\ASUS\\Desktop\\elktro\\outage_tool\\outage_qgis.py', encoding='utf-8').read())

PROFILE выбирает набор данных:
  'norge'  – вся электросеть Норвегии из FileGDB NVE (≈155 тыс. линий, 2090 подстанций, 734 тыс. опор);
  'region' – фрагмент сети (папка !!!ТОчняк, shp).
При первом запуске исходные слои приводятся к топологической сети по правилам (см. руководство)
и сохраняются в GeoPackage (слои net_edges, net_nodes, net_masts). Затем открывается панель.
Для всей Норвегии первое построение занимает несколько минут, последующие запуски – загрузка готовой сети.
"""
import math, time

# =========================== НАСТРОЙКИ ===========================
PROFILE = 'norge'          # 'norge' – вся сеть (FileGDB), 'region' – фрагмент (shp)
BASE = r'C:\Users\ASUS\Desktop\elktro'
PROFILES = {
    'norge': {'src': BASE + r'\Energi_0000_Norge_25833_Nettanlegg_FGDB\Energi_0000_Norge_25833_Nettanlegg_FGDB.gdb',
              'gpkg': BASE + r'\outage_tool\elektro_topo_norge.gpkg', 'tag': 'Норвегия'},
    'region': {'src': BASE + r'\!!!ТОчняк',
               'gpkg': BASE + r'\outage_tool\elektro_topo.gpkg', 'tag': 'район'},
}
BBOX = None        # (xmin, ymin, xmax, ymax) в EPSG:25833 – обработать только часть сети; None – всё
LOAD_MASTS = True  # False – не загружать опоры (быстрее и меньше памяти)
REBUILD = False    # True – перестроить топологию заново, даже если GeoPackage уже есть
# Параметры топологии (допуски в метрах) – словарь DEFAULT_PARAMS ниже
# =================================================================
SRC_DIR = PROFILES[PROFILE]['src']
OUT_GPKG = PROFILES[PROFILE]['gpkg']
TAG = PROFILES[PROFILE]['tag']

# ==== TOPO_CORE_BEGIN ====
# Чистый Python (без внешних библиотек): построение топологической сети
# и расчёт последствий аварийных отключений. Используется и в QGIS, и автономно.
from collections import defaultdict, deque

DEFAULT_PARAMS = {
    'SNAP_TOL': 10.0,        # м: концы линий ближе этого расстояния сливаются в один узел
    'TJUNC_TOL': 10.0,       # м: конец линии ближе этого к середине другой линии -> Т-ответвление
    'STATION_GROUP': 80.0,   # м: точки подстанции ближе этого объединяются в одну станцию
    'STATION_RADIUS': 100.0, # м: узлы линий в этом радиусе присоединяются к шинам станции
    'MAST_TOL': 25.0,        # м: опора привязывается к ближайшей линии в этом радиусе
    'MIN_LEN': 0.01,         # м: минимальная длина ребра
    'PARALLEL_COS': 0.94,    # |cos| угла > 0.94 (угол < 20°) - линии параллельны, Т-узел не создаётся
    'EXT_VOLT': 110.0,       # кВ: тупиковые концы линий >= этого напряжения = внешнее питание (граница района)
    'AUTO_FEED_RADIUS': 1500.0,  # м: изолированный участок РС подключается к ближайшей станции
                                 # условным кабельным вводом (0 - отключить правило)
}


def _d(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


def _plen(pts):
    return sum(_d(pts[i], pts[i + 1]) for i in range(len(pts) - 1))


class _UF:
    def __init__(self, n):
        self.p = list(range(n))

    def f(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def u(self, a, b):
        a, b = self.f(a), self.f(b)
        if a != b:
            self.p[max(a, b)] = min(a, b)


class _PointGrid:
    def __init__(self, pts, cell):
        self.cell = cell
        self.pts = pts
        self.g = defaultdict(list)
        for i, p in enumerate(pts):
            self.g[(int(p[0] // cell), int(p[1] // cell))].append(i)

    def near(self, p, r):
        c = self.cell
        cx0, cx1 = int((p[0] - r) // c), int((p[0] + r) // c)
        cy0, cy1 = int((p[1] - r) // c), int((p[1] + r) // c)
        out = []
        for cx in range(cx0, cx1 + 1):
            for cy in range(cy0, cy1 + 1):
                for i in self.g.get((cx, cy), ()):
                    if _d(self.pts[i], p) <= r:
                        out.append(i)
        return out


class _SegGrid:
    """Индекс отрезков полилиний: (edge_idx, seg_idx, a, b)."""

    def __init__(self, edges_pts, cell=50.0, skip=None):
        self.cell = cell
        self.g = defaultdict(list)
        self.segs = []
        for ei, pts in enumerate(edges_pts):
            if skip and skip(ei):
                continue
            for si in range(len(pts) - 1):
                a, b = pts[si], pts[si + 1]
                k = len(self.segs)
                self.segs.append((ei, si, a, b))
                # ячейки вдоль отрезка (шаг cell/2), а не весь охватывающий прямоугольник
                L = math.hypot(b[0] - a[0], b[1] - a[1])
                ns = int(L / (cell * 0.5)) + 1
                cells = set()
                for j in range(ns + 1):
                    t = j / ns
                    cells.add((int((a[0] + t * (b[0] - a[0])) // cell), int((a[1] + t * (b[1] - a[1])) // cell)))
                for c in cells:
                    self.g[c].append(k)

    def near(self, p, r):
        """Возвращает список (dist, edge_idx, seg_idx, t, proj) для отрезков в радиусе r."""
        c = self.cell
        seen = set()
        out = []
        rr = r + c * 0.3
        for cx in range(int((p[0] - rr) // c), int((p[0] + rr) // c) + 1):
            for cy in range(int((p[1] - rr) // c), int((p[1] + rr) // c) + 1):
                for k in self.g.get((cx, cy), ()):
                    if k in seen:
                        continue
                    seen.add(k)
                    ei, si, a, b = self.segs[k]
                    dx, dy = b[0] - a[0], b[1] - a[1]
                    L2 = dx * dx + dy * dy
                    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2))
                    q = (a[0] + t * dx, a[1] + t * dy)
                    dd = _d(p, q)
                    if dd <= r:
                        out.append((dd, ei, si, t, q))
        out.sort(key=lambda x: x[0])
        return out


def _volt_ok(v1, v2):
    return v1 is None or v2 is None or v1 == 0 or v2 == 0 or abs(v1 - v2) < 1e-6 or \
        {v1, v2} <= {22.0, 24.0} or {v1, v2} <= {11.0, 12.0}


def build_topology(lines, stations, masts, params=None, log=print):
    """
    lines:    [{'parts': [[(x,y),..],..], 'layer': str, 'src_id': int, 'name': str, 'volt': float|None,
                'level': str, 'owner': str}]
    stations: [{'pt': (x,y), 'name': str, 'volt': float|None, 'src_id': int}]
    masts:    [{'pt': (x,y), 'src_id': int}]
    Возвращает dict(edges=[...], nodes=[...], masts=[...], stats={...})
    """
    P = dict(DEFAULT_PARAMS)
    if params:
        P.update(params)
    st = {}
    t0 = time.perf_counter()

    # ---- Правило 1-2: одночастные геометрии, очистка, удаление дубликатов ----
    edges = []
    n_multi = n_short = n_dup = n_vtx = 0
    seen = set()
    for ln in lines:
        parts = ln['parts']
        if len(parts) > 1:
            n_multi += 1
        for part in parts:
            pts = []
            for p in part:
                p = (float(p[0]), float(p[1]))
                if pts and _d(pts[-1], p) < 1e-6:
                    n_vtx += 1
                    continue
                pts.append(p)
            if len(pts) < 2 or _plen(pts) < P['MIN_LEN']:
                n_short += 1
                continue
            key = tuple((round(x, 2), round(y, 2)) for x, y in pts)
            if key in seen or tuple(reversed(key)) in seen:
                n_dup += 1
                continue
            seen.add(key)
            e = dict(ln)
            e.pop('parts', None)
            e['pts'] = pts
            edges.append(e)
    st['source_lines'] = len(lines)
    st['multipart_exploded'] = n_multi
    st['short_removed'] = n_short
    st['duplicates_removed'] = n_dup
    st['dup_vertices_removed'] = n_vtx

    # ---- Правило 3: сведение концов линий в общие узлы (snapping) ----
    ends = []
    for e in edges:
        ends.append(e['pts'][0])
        ends.append(e['pts'][-1])
    grid = _PointGrid(ends, P['SNAP_TOL'])
    uf = _UF(len(ends))
    n_before = len(set(ends))
    for i, p in enumerate(ends):
        for j in grid.near(p, P['SNAP_TOL']):
            if j > i:
                uf.u(i, j)
    groups = defaultdict(list)
    for i in range(len(ends)):
        groups[uf.f(i)].append(i)
    node_xy = {}
    end_node = [None] * len(ends)
    nid = 0
    moved = 0
    for r, idx in groups.items():
        xs = sum(ends[i][0] for i in idx) / len(idx)
        ys = sum(ends[i][1] for i in idx) / len(idx)
        node_xy[nid] = (xs, ys)
        for i in idx:
            end_node[i] = nid
            if _d(ends[i], (xs, ys)) > 1e-6:
                moved += 1
        nid += 1
    for k, e in enumerate(edges):
        e['from'] = end_node[2 * k]
        e['to'] = end_node[2 * k + 1]
        e['pts'][0] = node_xy[e['from']]
        e['pts'][-1] = node_xy[e['to']]
    st['endpoints_snapped'] = moved
    st['nodes_after_snap'] = len(node_xy)

    # ---- Правило 4: Т-ответвления (конец линии на середине другой линии) ----
    inc = defaultdict(list)
    for k, e in enumerate(edges):
        inc[e['from']].append(k)
        inc[e['to']].append(k)
    sg = _SegGrid([e['pts'] for e in edges], cell=max(50.0, P['TJUNC_TOL'] * 3))
    splits = defaultdict(list)  # edge -> [(seg, t, node, proj)]
    for n, xy in node_xy.items():
        my = set(inc[n])
        if len(inc[n]) != 1:
            continue          # только «висячие» концы; узел-развилка не цепляется к проходящей линии
        k0 = inc[n][0]
        p0 = edges[k0]['pts']
        dv = (p0[1][0] - p0[0][0], p0[1][1] - p0[0][1]) if edges[k0]['from'] == n else \
             (p0[-2][0] - p0[-1][0], p0[-2][1] - p0[-1][1])
        myv = [edges[k].get('volt') for k in my]
        for dd, ei, si, t, q in sg.near(xy, P['TJUNC_TOL']):
            if ei in my:
                continue
            tgt = edges[ei]
            sa, sb = tgt['pts'][si], tgt['pts'][si + 1]
            sv = (sb[0] - sa[0], sb[1] - sa[1])
            nn = math.hypot(*dv) * math.hypot(*sv)
            if nn > 0 and abs(dv[0] * sv[0] + dv[1] * sv[1]) / nn > P['PARALLEL_COS']:
                continue      # почти параллельные линии (двухцепные ВЛ) не соединяем
            if _d(q, tgt['pts'][0]) < P['SNAP_TOL'] or _d(q, tgt['pts'][-1]) < P['SNAP_TOL']:
                continue
            if not any(_volt_ok(v, tgt.get('volt')) for v in myv):
                continue
            splits[ei].append((si, t, n, q))
            break
    # переносим узлы на проекцию и режем линии
    alias = {}
    for ei, lst in splits.items():
        for si, t, n, q in lst:
            node_xy[n] = q
    for k, e in enumerate(edges):
        e['pts'][0] = node_xy[e['from']]
        e['pts'][-1] = node_xy[e['to']]
    new_edges = []
    n_tj = 0
    for k, e in enumerate(edges):
        if k not in splits:
            new_edges.append(e)
            continue
        lst = sorted(splits[k], key=lambda s: (s[0], s[1]))
        pts = e['pts']
        cur = [pts[0]]
        cur_from = e['from']
        si_done = 0
        for si, t, n, q in lst:
            while si_done < si:
                si_done += 1
                cur.append(pts[si_done])
            if _d(cur[-1], q) > 1e-6:
                cur.append(q)
            if _plen(cur) < P['MIN_LEN']:
                alias[n] = cur_from
                cur = [q]
                continue
            ne = dict(e)
            ne['pts'] = cur
            ne['from'] = cur_from
            ne['to'] = n
            new_edges.append(ne)
            n_tj += 1
            cur = [q]
            cur_from = n
        while si_done < len(pts) - 1:
            si_done += 1
            cur.append(pts[si_done])
        ne = dict(e)
        ne['pts'] = cur
        ne['from'] = cur_from
        ne['to'] = e['to']
        if _plen(cur) >= P['MIN_LEN']:
            new_edges.append(ne)
        else:
            alias[e['to']] = cur_from
    for e in new_edges:
        e['from'] = alias.get(e['from'], e['from'])
        e['to'] = alias.get(e['to'], e['to'])
    edges = new_edges
    st['t_junctions'] = sum(len(v) for v in splits.values())

    # удаляем вырожденные петли
    loops = [e for e in edges if e['from'] == e['to'] and _plen(e['pts']) < 2 * P['SNAP_TOL']]
    edges = [e for e in edges if not (e['from'] == e['to'] and _plen(e['pts']) < 2 * P['SNAP_TOL'])]
    st['short_edges_collapsed'] = len(loops)
    for e in edges:
        e['kind'] = 'LINE'

    # ---- Правило 5: подстанции как узлы-источники с шинами ----
    spts = [s['pt'] for s in stations]
    suf = _UF(len(spts))
    sgrid = _PointGrid(spts, P['STATION_GROUP'])
    for i, p in enumerate(spts):
        for j in sgrid.near(p, P['STATION_GROUP']):
            if j > i:
                suf.u(i, j)
    sgroups = defaultdict(list)
    for i in range(len(spts)):
        sgroups[suf.f(i)].append(i)
    used_nodes = set()
    for e in edges:
        used_nodes.add(e['from'])
        used_nodes.add(e['to'])
    used_list = sorted(used_nodes)
    ngrid = _PointGrid([node_xy[n] for n in used_list], P['STATION_RADIUS'])
    station_nodes = {}
    best = {}  # line node -> (dist, station_node)
    for r, idx in sgroups.items():
        idx.sort(key=lambda i: (0 if 'TS' in (stations[i].get('name') or '').upper() else 1, i))
        main = stations[idx[0]]
        names = []
        for i in idx:
            nm = (stations[i].get('name') or '').strip()
            if nm and nm not in names:
                names.append(nm)
        volts = [stations[i].get('volt') for i in idx if stations[i].get('volt')]
        sn = nid
        nid += 1
        node_xy[sn] = main['pt']
        station_nodes[sn] = {'name': ' / '.join(names), 'volt': max(volts) if volts else None}
        for i in idx:
            for j in ngrid.near(spts[i], P['STATION_RADIUS']):
                ln = used_list[j]
                dd = _d(spts[i], node_xy[ln])
                if ln not in best or dd < best[ln][0]:
                    best[ln] = (dd, sn)
    n_bus = n_merge = 0
    merge = {}
    for ln, (dd, sn) in best.items():
        if dd < 0.5:
            merge[ln] = sn
            n_merge += 1
        else:
            edges.append({'kind': 'BUS', 'pts': [node_xy[sn], node_xy[ln]], 'from': sn, 'to': ln,
                          'layer': 'bus', 'name': 'Шина ' + station_nodes[sn]['name'],
                          'volt': None, 'level': '', 'owner': '', 'src_id': None})
            n_bus += 1
    for e in edges:
        e['from'] = merge.get(e['from'], e['from'])
        e['to'] = merge.get(e['to'], e['to'])
    st['stations'] = len(station_nodes)
    st['bus_edges'] = n_bus
    st['stations_isolated'] = sum(1 for sn in station_nodes if sn not in {b for _, b in best.values()})

    # ---- Правило 6: перенумерация, степени, типы узлов, внешнее питание ----
    deg = defaultdict(int)
    maxv = defaultdict(float)
    for e in edges:
        deg[e['from']] += 1
        deg[e['to']] += 1
        v = e.get('volt') or 0.0
        maxv[e['from']] = max(maxv[e['from']], v)
        maxv[e['to']] = max(maxv[e['to']], v)
    old2new = {}
    nodes = []
    for old in sorted(set(list(deg.keys()) + list(station_nodes.keys()))):
        new = len(nodes) + 1
        old2new[old] = new
        d = deg.get(old, 0)
        if old in station_nodes:
            ntype, name, volt, src = 'STATION', station_nodes[old]['name'], station_nodes[old]['volt'], 1
        else:
            ntype = 'END' if d == 1 else ('PASS' if d == 2 else 'BRANCH')
            name, volt, src = '', None, 0
            if d == 1 and maxv[old] >= P['EXT_VOLT'] > 0:
                ntype, name, volt, src = 'EXTERNAL', 'Внешнее питание (граница района)', maxv[old], 1
        nodes.append({'nid': new, 'x': node_xy[old][0], 'y': node_xy[old][1], 'ntype': ntype,
                      'degree': d, 'is_source': src, 'name': name, 'volt': volt})
    for k, e in enumerate(edges):
        e['eid'] = k + 1
        e['from'] = old2new[e['from']]
        e['to'] = old2new[e['to']]
        e['pts'][0] = (nodes[e['from'] - 1]['x'], nodes[e['from'] - 1]['y'])
        e['pts'][-1] = (nodes[e['to'] - 1]['x'], nodes[e['to'] - 1]['y'])
        e['length'] = _plen(e['pts'])

    # ---- Правило 7: связность ----
    def components():
        adj = defaultdict(list)
        for e in edges:
            adj[e['from']].append(e['to'])
            adj[e['to']].append(e['from'])
        comp = {}
        c = 0
        for n in nodes:
            if n['nid'] in comp:
                continue
            c += 1
            q = deque([n['nid']])
            comp[n['nid']] = c
            while q:
                u = q.popleft()
                for v in adj[u]:
                    if v not in comp:
                        comp[v] = c
                        q.append(v)
        return comp, c

    comp, c = components()
    fed = set(comp[n['nid']] for n in nodes if n['is_source'])
    st['components_initial'] = c
    st['components_with_source_initial'] = len(fed)
    # Правило 7б: условный кабельный ввод к изолированным участкам распределительной сети
    n_virt = 0
    if P['AUTO_FEED_RADIUS'] > 0:
        stn = [n for n in nodes if n['ntype'] == 'STATION']
        bycomp = defaultdict(list)
        for n in nodes:
            if comp[n['nid']] not in fed and n['ntype'] != 'STATION':
                bycomp[comp[n['nid']]].append(n)
        sgr = _PointGrid([(s_['x'], s_['y']) for s_ in stn], max(P['AUTO_FEED_RADIUS'], 1.0))
        for cc, ns in bycomp.items():
            bestv = None
            for a in ns:
                for j in sgr.near((a['x'], a['y']), P['AUTO_FEED_RADIUS']):
                    s_ = stn[j]
                    dd = math.hypot(a['x'] - s_['x'], a['y'] - s_['y'])
                    if bestv is None or dd < bestv[0]:
                        bestv = (dd, a, s_)
            if bestv:
                dd, a, s_ = bestv
                edges.append({'kind': 'VIRT', 'pts': [(s_['x'], s_['y']), (a['x'], a['y'])],
                              'from': s_['nid'], 'to': a['nid'], 'layer': 'virtual',
                              'name': 'Условный ввод от ' + s_['name'], 'volt': None, 'level': '',
                              'owner': '', 'src_id': None, 'eid': len(edges) + 1, 'length': dd})
                a['degree'] += 1
                s_['degree'] += 1
                n_virt += 1
        comp, c = components()
        fed = set(comp[n['nid']] for n in nodes if n['is_source'])
    st['virtual_feeds'] = n_virt
    for n in nodes:
        n['comp'] = comp[n['nid']]
        n['fed'] = 1 if n['comp'] in fed else 0
    for e in edges:
        e['comp'] = comp[e['from']]
        e['fed'] = 1 if e['comp'] in fed else 0
    st['components'] = c
    st['components_with_source'] = len(fed)
    st['edges'] = len(edges)
    st['line_edges'] = sum(1 for e in edges if e['kind'] == 'LINE')
    st['nodes'] = len(nodes)
    st['node_types'] = dict((t, sum(1 for n in nodes if n['ntype'] == t)) for t in ('STATION', 'EXTERNAL', 'BRANCH', 'PASS', 'END'))
    st['len_km'] = sum(e['length'] for e in edges if e['kind'] == 'LINE') / 1000.0
    st['len_fed_km'] = sum(e['length'] for e in edges if e['kind'] == 'LINE' and e['fed']) / 1000.0

    # ---- Правило 8: привязка опор к линиям ----
    lsg = _SegGrid([e['pts'] for e in edges], cell=50.0, skip=lambda i: edges[i]['kind'] != 'LINE')
    out_m = []
    n_bound = 0
    for k, m in enumerate(masts):
        near = lsg.near(m['pt'], P['MAST_TOL'])
        if near:
            dd, ei = near[0][0], near[0][1]
            eid = edges[ei]['eid']
            n_bound += 1
        else:
            dd, eid = None, None
        out_m.append({'mid': k + 1, 'pt': m['pt'], 'edge_id': eid, 'dist': dd, 'src_id': m.get('src_id')})
    st['masts'] = len(masts)
    st['masts_bound'] = n_bound
    st['build_sec'] = time.perf_counter() - t0
    return {'edges': edges, 'nodes': nodes, 'masts': out_m, 'stats': st, 'params': P}


class OutageGraph:
    """Граф сети для расчёта отключений. edges: [(eid, from, to)], sources: [nid]."""
    ST_ON, ST_OFF, ST_DEAD, ST_NOFEED = 0, 1, 2, 3

    def __init__(self, edges, node_ids, sources):
        self.edges = {eid: (a, b) for eid, a, b in edges}
        self.node_ids = set(node_ids)
        self.sources = set(sources)
        self.adj = defaultdict(list)
        for eid, a, b in edges:
            self.adj[a].append((eid, b))
            self.adj[b].append((eid, a))
        self.base_nodes = self._reach(set(), set())
        self._base_es = self._base_ns = None

    def _reach(self, off_e, off_n):
        seen = set()
        q = deque()
        for s in self.sources:
            if s not in off_n:
                seen.add(s)
                q.append(s)
        adj = self.adj
        while q:
            u = q.popleft()
            for eid, v in adj[u]:
                if eid in off_e or v in off_n or v in seen:
                    continue
                seen.add(v)
                q.append(v)
        return seen

    def analyze(self, off_edges=(), off_nodes=()):
        """Возвращает (edge_state, node_state, secs). Состояния:
        0 - под напряжением, 1 - отключено оператором, 2 - обесточено в результате отключения,
        3 - не имеет связи с источником даже в нормальном режиме.
        Инкрементальный алгоритм: пересчитываются только фрагменты сети, примыкающие к отключённым
        элементам; поиск в ширину из фрагмента останавливается при достижении включённого источника."""
        t0 = time.perf_counter()
        off_e, off_n = set(off_edges), set(off_nodes)
        base = self.base_nodes
        if self._base_es is None:
            self._base_es = {eid: (0 if (a in base or b in base) else 3) for eid, (a, b) in self.edges.items()}
            self._base_ns = {n: (0 if n in base else 3) for n in self.node_ids}
        es = dict(self._base_es)
        ns = dict(self._base_ns)
        self.touched_e, self.touched_n = set(), set()   # элементы, состояние которых может отличаться от нормального
        if not off_e and not off_n:
            return es, ns, time.perf_counter() - t0
        adj, src = self.adj, self.sources
        # стартовые узлы: концы отключённых рёбер и соседи отключённых узлов
        start = set()
        for eid in off_e:
            a, b = self.edges[eid]
            start.add(a); start.add(b)
        for n in off_n:
            for eid, v in adj[n]:
                start.add(v)
        start -= off_n
        live, dead = set(), set()
        for s0 in start:
            if s0 in live or s0 in dead or s0 not in base:
                continue
            seen = {s0}
            q = deque([s0])
            found = s0 in src
            while q and not found:
                u = q.popleft()
                for eid, v in adj[u]:
                    if eid in off_e or v in off_n or v in seen:
                        continue
                    if v in live or v in src:
                        found = True
                        break
                    seen.add(v)
                    q.append(v)
            if found:
                live |= seen
            else:
                dead |= seen
        for n in off_n:
            ns[n] = 1
        for n in dead:
            ns[n] = 2
        touched = set(off_e)
        for n in dead | off_n:
            for eid, v in adj[n]:
                touched.add(eid)
        self.touched_e, self.touched_n = touched, dead | off_n
        for eid in touched:
            if eid in off_e:
                es[eid] = 1
                continue
            a, b = self.edges[eid]
            la = a in base and a not in dead and a not in off_n
            lb = b in base and b not in dead and b not in off_n
            if la or lb:
                es[eid] = 0
            elif a in base or b in base:
                es[eid] = 2
            else:
                es[eid] = 3
        return es, ns, time.perf_counter() - t0

    def analyze_full(self, off_edges=(), off_nodes=()):
        """Эталонный (полный) расчёт многоисточниковым поиском в ширину – для проверки."""
        t0 = time.perf_counter()
        off_e, off_n = set(off_edges), set(off_nodes)
        live = self._reach(off_e, off_n)
        base = self.base_nodes
        es = {}
        for eid, (a, b) in self.edges.items():
            if eid in off_e:
                es[eid] = 1
            elif a in live or b in live:
                es[eid] = 0
            elif a in base or b in base:
                es[eid] = 2
            else:
                es[eid] = 3
        ns = {}
        for n in self.node_ids:
            if n in off_n:
                ns[n] = 1
            elif n in live:
                ns[n] = 0
            elif n in base:
                ns[n] = 2
            else:
                ns[n] = 3
        return es, ns, time.perf_counter() - t0
# ==== TOPO_CORE_END ====

# =====================================================================
#                ЧАСТЬ QGIS: построение слоёв, стили, панель
# =====================================================================
import os, random, csv, datetime
from collections import defaultdict
from qgis.core import (QgsProject, QgsVectorLayer, QgsFeature, QgsGeometry, QgsPointXY, QgsVectorFileWriter,
                       QgsCoordinateTransform, QgsRuleBasedRenderer, QgsLineSymbol, QgsMarkerSymbol,
                       QgsPalLayerSettings, QgsTextFormat, QgsVectorLayerSimpleLabeling, QgsRectangle,
                       QgsFeatureRequest, NULL)
from qgis.PyQt.QtCore import Qt
from qgis.PyQt.QtGui import QColor, QFont
from qgis.PyQt.QtWidgets import (QDockWidget, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QListWidget,
                                 QLabel, QTextEdit, QCheckBox, QMessageBox, QGroupBox, QLineEdit, QApplication)
from qgis.gui import QgsMapTool
from qgis.utils import iface

SRC_LAYERS = ['el_luftlinje', 'el_hengekabel', 'el_sjokabel']
STATION_LAYER = 'el_transformatorstasjon'
MAST_LAYER = 'el_mast'
L_EDGES, L_NODES, L_MASTS = 'net_edges', 'net_nodes', 'net_masts'
TITLE = {L_EDGES: 'Сеть (%s): линии (рёбра)' % TAG, L_NODES: 'Сеть (%s): узлы' % TAG, L_MASTS: 'Сеть (%s): опоры' % TAG}
GROUP_NAME = 'Анализ аварийных отключений – %s' % TAG
STATE_TXT = {0: 'под напряжением', 1: 'отключено оператором', 2: 'обесточено', 3: 'нет связи с источником'}
PICK_MAX_TOL = 1500.0   # м: при более мелком масштабе выбор на карте отключается (слишком много объектов)


def _msg(text):
    print(text)
    try:
        iface.mainWindow().statusBar().showMessage(text)
        QApplication.processEvents()
    except Exception:
        pass


def _val(f, name, conv=None):
    try:
        v = f.attribute(name)
    except KeyError:
        return None
    if v is None or v == NULL:
        return None
    if conv:
        try:
            return conv(v)
        except (TypeError, ValueError):
            return None
    return v


def _src_layer(name, src):
    """Слой исходных данных: из FileGDB (.gdb), из папки shp или уже открытый в проекте из того же источника."""
    s_norm = src.replace('\\', '/').lower()
    for l in QgsProject.instance().mapLayersByName(name):
        if isinstance(l, QgsVectorLayer) and l.isValid() and l.source().replace('\\', '/').lower().startswith(s_norm):
            return l
    if src.lower().endswith('.gdb'):
        l = QgsVectorLayer('%s|layername=%s' % (src, name), name, 'ogr')
        return l if l.isValid() else None
    p = os.path.join(src, name + '.shp')
    if os.path.exists(p):
        l = QgsVectorLayer(p, name, 'ogr')
        if l.isValid():
            return l
    return None


def _request(lyr, fields, bbox=None, geom=True):
    rq = QgsFeatureRequest()
    names = [n for n in fields if lyr.fields().indexOf(n) >= 0]
    rq.setSubsetOfAttributes(names, lyr.fields())
    if not geom:
        rq.setFlags(QgsFeatureRequest.NoGeometry)
    if bbox:
        rq.setFilterRect(QgsRectangle(*bbox))
    return rq


def _xy_parts(g):
    if g is None or g.isEmpty():
        return []
    if g.isMultipart():
        return [[(p.x(), p.y()) for p in part] for part in g.asMultiPolyline()]
    return [[(p.x(), p.y()) for p in g.asPolyline()]]


def _xy_point(g):
    if g is None or g.isEmpty():
        return None
    p = g.asMultiPoint()[0] if g.isMultipart() else g.asPoint()
    return (p.x(), p.y())


def _write_layer(mem, path, lname, first):
    o = QgsVectorFileWriter.SaveVectorOptions()
    o.driverName = 'GPKG'
    o.layerName = lname
    o.fileEncoding = 'UTF-8'
    o.actionOnExistingFile = (QgsVectorFileWriter.CreateOrOverwriteFile if first
                              else QgsVectorFileWriter.CreateOrOverwriteLayer)
    res = QgsVectorFileWriter.writeAsVectorFormatV3(mem, path, QgsProject.instance().transformContext(), o)
    if res[0] != QgsVectorFileWriter.NoError:
        raise RuntimeError('Ошибка записи %s: %s' % (lname, res[1]))


def _add_chunked(layer, feats, chunk=50000):
    for i in range(0, len(feats), chunk):
        layer.dataProvider().addFeatures(feats[i:i + chunk])


def build_network(src=None, gpkg=None, params=None, load=True, bbox=None, masts=None):
    """Правила топологии -> GeoPackage со слоями net_edges / net_nodes / net_masts."""
    src = src or SRC_DIR
    gpkg = gpkg or OUT_GPKG
    bbox = bbox if bbox is not None else BBOX
    masts_on = LOAD_MASTS if masts is None else masts
    os.makedirs(os.path.dirname(gpkg), exist_ok=True)
    t0 = time.perf_counter()
    lines, crs = [], None
    for name in SRC_LAYERS:
        lyr = _src_layer(name, src)
        if lyr is None:
            print('  ! слой не найден:', name)
            continue
        crs = crs or lyr.crs()
        _msg('Чтение %s…' % name)
        for f in lyr.getFeatures(_request(lyr, ['OBJECTID', 'navn', 'spenning', 'nettniva', 'eier'], bbox)):
            parts = _xy_parts(f.geometry())
            if parts:
                lines.append({'parts': parts, 'layer': name, 'src_id': _val(f, 'OBJECTID', int) or f.id(),
                              'name': _val(f, 'navn', str) or '', 'volt': _val(f, 'spenning', float),
                              'level': _val(f, 'nettniva', str) or '', 'owner': _val(f, 'eier', str) or ''})
    if not lines:
        raise RuntimeError('Не найдены исходные линии. Проверьте источник: %s' % src)
    stations, masts = [], []
    lyr = _src_layer(STATION_LAYER, src)
    if lyr:
        for f in lyr.getFeatures(_request(lyr, ['OBJECTID', 'navn', 'spenning'], bbox)):
            p = _xy_point(f.geometry())
            if p:
                v = _val(f, 'spenning', float)
                stations.append({'pt': p, 'name': _val(f, 'navn', str) or '', 'volt': v or None,
                                 'src_id': _val(f, 'OBJECTID', int) or f.id()})
    lyr = _src_layer(MAST_LAYER, src) if masts_on else None
    if lyr:
        _msg('Чтение опор…')
        for f in lyr.getFeatures(_request(lyr, ['OBJECTID'], bbox)):
            p = _xy_point(f.geometry())
            if p:
                masts.append({'pt': p, 'src_id': _val(f, 'OBJECTID', int) or f.id()})
    _msg('Прочитано за %.0f с: линий %d, станций %d, опор %d. Построение топологии…'
         % (time.perf_counter() - t0, len(lines), len(stations), len(masts)))
    res = build_topology(lines, stations, masts, params)
    del lines
    authid = crs.authid() if crs else 'EPSG:25833'
    _msg('Топология построена (%.0f с). Запись GeoPackage…' % (time.perf_counter() - t0))

    me = QgsVectorLayer('LineString?crs=%s&field=edge_id:integer&field=from_node:integer&field=to_node:integer'
                        '&field=kind:string(8)&field=src_layer:string(24)&field=src_id:integer&field=name:string(120)'
                        '&field=volt:double&field=level:string(40)&field=owner:string(80)&field=length:double'
                        '&field=comp:integer&field=fed:integer&field=state:integer' % authid, L_EDGES, 'memory')
    feats = []
    for e in res['edges']:
        f = QgsFeature(me.fields())
        f.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(x, y) for x, y in e['pts']]))
        f.setAttributes([e['eid'], e['from'], e['to'], e['kind'], e.get('layer'), e.get('src_id'),
                         (e.get('name') or '')[:120], e.get('volt'), e.get('level'), (e.get('owner') or '')[:80],
                         round(e['length'], 2), e['comp'], e['fed'], 0 if e['fed'] else 3])
        feats.append(f)
    _add_chunked(me, feats)

    mn = QgsVectorLayer('Point?crs=%s&field=node_id:integer&field=ntype:string(10)&field=degree:integer'
                        '&field=is_source:integer&field=name:string(120)&field=volt:double&field=comp:integer'
                        '&field=fed:integer&field=state:integer' % authid, L_NODES, 'memory')
    feats = []
    for n in res['nodes']:
        f = QgsFeature(mn.fields())
        f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(n['x'], n['y'])))
        f.setAttributes([n['nid'], n['ntype'], n['degree'], n['is_source'], (n['name'] or '')[:120], n['volt'],
                         n['comp'], n['fed'], 0 if n['fed'] else 3])
        feats.append(f)
    _add_chunked(mn, feats)

    fed_e = {e['eid']: e['fed'] for e in res['edges']}
    mm = QgsVectorLayer('Point?crs=%s&field=mast_id:integer&field=src_id:integer&field=edge_id:integer'
                        '&field=dist:double&field=state:integer' % authid, L_MASTS, 'memory')
    feats = []
    for m in res['masts']:
        f = QgsFeature(mm.fields())
        f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(*m['pt'])))
        st = 0 if (m['edge_id'] and fed_e.get(m['edge_id'])) else 3
        f.setAttributes([m['mid'], m['src_id'], m['edge_id'],
                         round(m['dist'], 2) if m['dist'] is not None else None, st])
        feats.append(f)
    _add_chunked(mm, feats)
    feats = None

    # снимаем слои с карты, если они уже открыты (иначе файл заблокирован)
    for lid, l in list(QgsProject.instance().mapLayers().items()):
        if l.source().replace('\\', '/').lower().startswith(gpkg.replace('\\', '/').lower()):
            QgsProject.instance().removeMapLayer(lid)
    _write_layer(me, gpkg, L_EDGES, True)
    _write_layer(mn, gpkg, L_NODES, False)
    _write_layer(mm, gpkg, L_MASTS, False)
    me = mn = mm = None

    st = res['stats']
    rep = os.path.splitext(gpkg)[0] + '_topology_report.txt'
    with open(rep, 'w', encoding='utf-8') as fh:
        fh.write('Отчёт о построении топологии (%s), источник: %s\n' % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), src))
        fh.write('  BBOX = %s\n' % (bbox,))
        for k, v in res['params'].items():
            fh.write('  параметр %s = %s\n' % (k, v))
        for k, v in st.items():
            fh.write('  %s: %s\n' % (k, v))
        fh.write('  total_sec (с чтением и записью): %.1f\n' % (time.perf_counter() - t0))
    _msg('Топология построена за %.0f с: рёбер %d, узлов %d, компонент %d (с источником %d). Файл: %s'
         % (time.perf_counter() - t0, st['edges'], st['nodes'], st['components'], st['components_with_source'], gpkg))
    if load:
        load_network_layers(gpkg)
    return st


# ------------------------------------------------------------------ стили
def _line(color, width, style='solid'):
    return QgsLineSymbol.createSimple({'color': color, 'width': str(width), 'line_style': style,
                                       'capstyle': 'round', 'joinstyle': 'round'})


def _mark(shape, color, size, outline='#222222'):
    return QgsMarkerSymbol.createSimple({'name': shape, 'color': color, 'size': str(size),
                                         'outline_color': outline, 'outline_width': '0.2'})


def _rules(layer, rules):
    root = QgsRuleBasedRenderer.Rule(None)
    for sym, expr, label, min_scale in rules:
        if min_scale == -1:     # правило только для мелких масштабов (мельче 1:250 000)
            r = QgsRuleBasedRenderer.Rule(sym, 250000, 0, expr, label)
        else:
            r = QgsRuleBasedRenderer.Rule(sym, 0, min_scale, expr, label)
        root.appendChild(r)
    layer.setRenderer(QgsRuleBasedRenderer(root))


def style_layers(le, ln, lm):
    if le:
        _rules(le, [
            (_line('#bdbdbd', 0.26), '"state" = 3', 'Нет связи с источником (разрыв в исходных данных)', 0),
            (_line('#666666', 0.3, 'dot'), '"state" = 0 AND "kind" <> \'LINE\'', 'Шина / условный ввод', 0),
            (_line('#2e8b57', 0.35), '"state" = 0 AND "kind" = \'LINE\' AND coalesce("volt",0) < 50', 'до 50 кВ — в работе', 0),
            (_line('#8e44ad', 0.6), '"state" = 0 AND "kind" = \'LINE\' AND "volt" >= 50 AND "volt" < 110', '50–109 кВ — в работе', 0),
            (_line('#1f3a93', 0.8), '"state" = 0 AND "kind" = \'LINE\' AND "volt" >= 110', '110–420 кВ — в работе', 0),
            (_line('#ff8c00', 1.0), '"state" = 2', 'ОБЕСТОЧЕНО вследствие отключения', 0),
            (_line('#e00000', 1.5, 'dash'), '"state" = 1', 'ОТКЛЮЧЕНО оператором', 0),
        ])
    if lm:
        _rules(lm, [
            (_mark('circle', '#9e9e9e', 0.7, '#9e9e9e'), '"state" IN (0, 3)', 'Опора', 60000),
            (_mark('circle', '#b35900', 0.8, '#b35900'), '"state" = 2', 'Опора обесточенной линии', 60000),
            (_mark('circle', '#e00000', 0.9, '#e00000'), '"state" = 1', 'Опора отключённой линии', 60000),
        ])
    if ln:
        _rules(ln, [
            (_mark('circle', '#424242', 1.0), '"state" IN (0, 3) AND "is_source" = 0', 'Узел сети', 60000),
            (_mark('circle', '#ff8c00', 1.2, '#ff8c00'), '"state" = 2', 'Обесточенный узел', 250000),
            (_mark('triangle', '#6a5acd', 3.0), '"state" = 0 AND "is_source" = 1 AND "ntype" <> \'STATION\'',
             'Внешнее питание / электростанция (тупик ВЛ ≥110 кВ)', 250000),
            (_mark('square', '#0050c8', 3.6), '"state" = 0 AND "ntype" = \'STATION\'', 'Трансформаторная станция (источник)', 250000),
            (_mark('square', '#0050c8', 0.7, '#0050c8'), '"state" = 0 AND "ntype" = \'STATION\'', 'Станция (мелкий масштаб)', -1),
            (_mark('cross2', '#e00000', 5.0, '#e00000'), '"state" = 1', 'ОТКЛЮЧЁННЫЙ узел', 0),
        ])
        s = QgsPalLayerSettings()
        s.fieldName = 'CASE WHEN "ntype" = \'STATION\' THEN "name" END'
        s.isExpression = True
        s.scaleVisibility = True
        s.minimumScale = 400000   # подписи станций – только при масштабе крупнее 1:400 000
        s.maximumScale = 0
        tf = QgsTextFormat()
        tf.setSize(8)
        tf.setColor(QColor('#0a2a6b'))
        s.setFormat(tf)
        ln.setLabeling(QgsVectorLayerSimpleLabeling(s))
        ln.setLabelsEnabled(True)
    for l in (le, ln, lm):
        if l:
            l.triggerRepaint()


def load_network_layers(gpkg=None):
    gpkg = gpkg or OUT_GPKG
    prj = QgsProject.instance()
    root = prj.layerTreeRoot()
    grp = root.findGroup(GROUP_NAME) or root.insertGroup(0, GROUP_NAME)
    out = {}
    for lname in (L_NODES, L_EDGES, L_MASTS):
        for l in prj.mapLayersByName(TITLE[lname]):
            prj.removeMapLayer(l.id())
        l = QgsVectorLayer('%s|layername=%s' % (gpkg, lname), TITLE[lname], 'ogr')
        if not l.isValid():
            raise RuntimeError('Не удалось открыть %s из %s' % (lname, gpkg))
        prj.addMapLayer(l, False)
        grp.addLayer(l)
        out[lname] = l
    style_layers(out[L_EDGES], out[L_NODES], out[L_MASTS])
    return out


# ------------------------------------------------------------------ инструмент карты
class PickTool(QgsMapTool):
    def __init__(self, canvas, callback):
        super().__init__(canvas)
        self.callback = callback
        self.setCursor(Qt.CrossCursor)

    def canvasReleaseEvent(self, e):
        pt = self.toMapCoordinates(e.pos())
        tol = self.canvas().mapUnitsPerPixel() * 7
        self.callback(pt, tol, e.button() == Qt.RightButton)


# ------------------------------------------------------------------ анализатор
class OutageTool:
    def __init__(self):
        prj = QgsProject.instance()

        def get(lname):
            ls = prj.mapLayersByName(TITLE[lname])
            return ls[0] if ls else None
        self.le, self.ln, self.lm = get(L_EDGES), get(L_NODES), get(L_MASTS)
        if not self.le or not self.ln:
            raise RuntimeError('Слои сети не загружены. Выполните build_network() или load_network_layers().')
        self.off_e, self.off_n = set(), set()
        self.last = None
        self.manual_t0 = None
        self._node_idx = self._edge_idx = None
        self.reload()

    # ---- загрузка графа
    def reload(self):
        t0 = time.perf_counter()
        _msg('Загрузка графа сети…')
        self.E, self.N = {}, {}
        self.cur_e, self.cur_n = {}, {}
        if self._reload_sqlite():
            pass
        else:
            self._reload_qgis()
        self.fid2eid = {v['fid']: k for k, v in self.E.items()}
        self.fid2nid = {v['fid']: k for k, v in self.N.items()}
        self._reload_masts()
        self.G = OutageGraph([(k, v['a'], v['b']) for k, v in self.E.items()], list(self.N.keys()),
                             [k for k, v in self.N.items() if v['src']])
        _msg('Граф загружен: %d рёбер, %d узлов, %d источников, %d опор, %.1f с'
             % (len(self.E), len(self.N), len(self.G.sources), self.n_masts, time.perf_counter() - t0))

    def _reload_sqlite(self):
        """Быстрое чтение атрибутов напрямую из GeoPackage (SQLite)."""
        import sqlite3
        src_e, src_n = self.le.source().split('|')[0], self.ln.source().split('|')[0]
        if not (src_e.lower().endswith('.gpkg') and src_n.lower().endswith('.gpkg')):
            return False
        try:
            con = sqlite3.connect(src_e)
            for fid, eid, a, b, kind, name, volt, lay, ln, st in con.execute(
                    'SELECT fid, edge_id, from_node, to_node, kind, name, volt, src_layer, length, state FROM %s' % L_EDGES):
                self.E[eid] = {'fid': fid, 'a': a, 'b': b, 'kind': kind, 'name': name or '', 'volt': volt,
                               'layer': lay or '', 'len': ln or 0.0}
                self.cur_e[eid] = st
            con.close()
            con = sqlite3.connect(src_n)
            for fid, nid, src, nt, name, st in con.execute(
                    'SELECT fid, node_id, is_source, ntype, name, state FROM %s' % L_NODES):
                self.N[nid] = {'fid': fid, 'src': src or 0, 'type': nt, 'name': name or ''}
                self.cur_n[nid] = st
            con.close()
            return True
        except Exception as ex:
            print('Чтение через SQLite невозможно (%s), используется QGIS' % ex)
            self.E, self.N, self.cur_e, self.cur_n = {}, {}, {}, {}
            return False

    def _reload_qgis(self):
        for f in self.le.getFeatures(_request(self.le, ['edge_id', 'from_node', 'to_node', 'kind', 'name', 'volt',
                                                        'src_layer', 'length', 'state'], geom=False)):
            eid = f['edge_id']
            self.E[eid] = {'fid': f.id(), 'a': f['from_node'], 'b': f['to_node'], 'kind': f['kind'],
                           'name': _val(f, 'name', str) or '', 'volt': _val(f, 'volt', float),
                           'layer': _val(f, 'src_layer', str) or '', 'len': _val(f, 'length', float) or 0.0}
            self.cur_e[eid] = _val(f, 'state', int)
        for f in self.ln.getFeatures(_request(self.ln, ['node_id', 'is_source', 'ntype', 'name', 'state'], geom=False)):
            nid = f['node_id']
            self.N[nid] = {'fid': f.id(), 'src': _val(f, 'is_source', int) or 0, 'type': f['ntype'],
                           'name': _val(f, 'name', str) or ''}
            self.cur_n[nid] = _val(f, 'state', int)

    def _reload_masts(self):
        # опоры не загружаются в память: их идентификаторы выбираются из GeoPackage по индексу edge_id
        self.mdb = None
        self.n_masts = 0
        if self.lm:
            try:
                import sqlite3
                path = self.lm.source().split('|')[0]
                self.mdb = sqlite3.connect(path, check_same_thread=False)
                self.mdb.execute('CREATE INDEX IF NOT EXISTS idx_%s_edge ON %s(edge_id)' % (L_MASTS, L_MASTS))
                self.mdb.commit()
                self.n_masts = self.mdb.execute('SELECT count(*) FROM %s' % L_MASTS).fetchone()[0]
            except Exception as ex:
                print('Опоры: доступ через SQLite невозможен (%s)' % ex)
                self.mdb = None

    def _masts_of(self, eids, count=False):
        """Идентификаторы (или число) опор, привязанных к рёбрам eids."""
        if not self.mdb or not eids:
            return 0 if count else []
        eids = list(eids)
        res = 0 if count else []
        for i in range(0, len(eids), 5000):
            lst = ','.join(str(int(e)) for e in eids[i:i + 5000])
            if count:
                res += self.mdb.execute('SELECT count(*) FROM %s WHERE edge_id IN (%s)' % (L_MASTS, lst)).fetchone()[0]
            else:
                res += [(r[0], r[1]) for r in self.mdb.execute('SELECT fid, edge_id FROM %s WHERE edge_id IN (%s)' % (L_MASTS, lst))]
        return res

    # ---- управление списком отключений
    def disconnect(self, lines=(), nodes=(), run=True):
        """Добавить линии (edge_id) и узлы (node_id) к отключению и пересчитать. Возвращает краткий итог."""
        bad = [e for e in lines if e not in self.E] + [n for n in nodes if n not in self.N]
        self.off_e.update(e for e in lines if e in self.E)
        self.off_n.update(n for n in nodes if n in self.N)
        if bad:
            print('Не найдены id:', bad)
        if not run:
            return None
        r = self.run()
        dead = [k for k, s in r['es'].items() if s == 2 and self.E[k]['kind'] == 'LINE']
        km = sum(self.E[k]['len'] for k in dead) / 1000
        return 'Обесточено участков: %d (%.2f км), расчёт %.1f мс' % (len(dead), km, r['t_calc'] * 1000)

    def toggle_line(self, eid):
        self.off_e.symmetric_difference_update({eid})

    def toggle_node(self, nid):
        self.off_n.symmetric_difference_update({nid})

    def clear(self):
        self.off_e.clear()
        self.off_n.clear()

    def add_selected(self):
        n = 0
        for f in self.le.selectedFeatures():
            self.off_e.add(f['edge_id'])
            n += 1
        for f in self.ln.selectedFeatures():
            self.off_n.add(f['node_id'])
            n += 1
        self.le.removeSelection()
        self.ln.removeSelection()
        return n

    # ---- расчёт
    def run(self, paint=True):
        es, ns, dt = self.G.analyze(self.off_e, self.off_n)
        t1 = time.perf_counter()
        if paint:
            self._paint(es, ns)
        self.last = {'es': es, 'ns': ns, 't_calc': dt, 't_total': dt + (time.perf_counter() - t1)}
        return self.last

    def _paint(self, es, ns):
        # после первого расчёта сравниваются только элементы, затронутые текущим и предыдущим отключением
        te, tn = getattr(self.G, 'touched_e', None), getattr(self.G, 'touched_n', None)
        if te is not None and getattr(self, '_prev_te', None) is not None:
            keys_e, keys_n = te | self._prev_te, tn | self._prev_tn
        else:
            keys_e, keys_n = es.keys(), ns.keys()
        self._prev_te, self._prev_tn = set(te or ()), set(tn or ())
        ch_edges = [k for k in keys_e if self.cur_e.get(k) != es[k]]
        ms = {fid: es[eid] for fid, eid in self._masts_of(ch_edges)}
        for lyr, cur, new, keys, key in ((self.le, self.cur_e, es, keys_e, lambda k: self.E[k]['fid']),
                                         (self.ln, self.cur_n, ns, keys_n, lambda k: self.N[k]['fid'])):
            idx = lyr.fields().indexOf('state')
            ch = {}
            for k in keys:
                s = new[k]
                if cur.get(k) != s:
                    ch[key(k)] = {idx: s}
                    cur[k] = s
            if ch:
                lyr.dataProvider().changeAttributeValues(ch)
                lyr.triggerRepaint()
        if self.lm and ms:
            idx = self.lm.fields().indexOf('state')
            self.lm.dataProvider().changeAttributeValues({fid: {idx: s} for fid, s in ms.items()})
            self.lm.triggerRepaint()

    def reset(self):
        self.clear()
        return self.run()

    # ---- отчёт
    def report(self):
        if not self.last:
            self.run()
        es, ns = self.last['es'], self.last['ns']
        dead = [k for k, s in es.items() if s == 2 and self.E[k]['kind'] == 'LINE']
        by_v = {}
        for k in dead:
            v = self.E[k]['volt']
            key = '%g кВ' % v if v else 'н/д'
            by_v[key] = by_v.get(key, 0.0) + self.E[k]['len']
        names = sorted({self.E[k]['name'] for k in dead if self.E[k]['name']})
        dnodes = [k for k, s in ns.items() if s == 2]
        dst = [self.N[k]['name'] for k in dnodes if self.N[k]['type'] == 'STATION']
        n_masts = self._masts_of([k for k, s in es.items() if s == 2], count=True)
        lines = ['<b>Отключено оператором:</b> линий %d, узлов %d' % (len(self.off_e), len(self.off_n)),
                 '<b style="color:#d35400">Обесточено участков линий:</b> %d (%.2f км)'
                 % (len(dead), sum(self.E[k]['len'] for k in dead) / 1000.0)]
        for k in sorted(by_v, key=lambda x: -by_v[x]):
            lines.append('&nbsp;&nbsp;• %s: %.2f км' % (k, by_v[k] / 1000.0))
        lines.append('<b>Обесточено узлов:</b> %d (ответвлений %d, конечных %d)'
                     % (len(dnodes), sum(1 for k in dnodes if self.N[k]['type'] == 'BRANCH'),
                        sum(1 for k in dnodes if self.N[k]['type'] == 'END')))
        lines.append('<b>Обесточено опор:</b> %d' % n_masts)
        if dst:
            lines.append('<b>Потеряли питание подстанции (%d):</b> %s' % (len(dst), '; '.join(dst[:20])))
        off_src = [self.N[k]['name'] or self.N[k]['type'] for k in self.off_n if self.N[k]['src']]
        if off_src:
            lines.append('<b>Отключены источники:</b> ' + '; '.join(off_src))
        if names:
            lines.append('<b>Затронутые линии (%d):</b> %s%s' % (len(names), '; '.join(names[:25]),
                                                                 ' …' if len(names) > 25 else ''))
        lines.append('<i>Время расчёта: %.1f мс (с перерисовкой %.1f мс); сеть: %d рёбер</i>'
                     % (self.last['t_calc'] * 1000, self.last['t_total'] * 1000, len(self.E)))
        return '<br>'.join(lines)

    def export_csv(self, path=None):
        if not self.last:
            self.run()
        path = path or os.path.join(os.path.dirname(OUT_GPKG), 'reports',
                                    'outage_%s_%s.csv' % (PROFILE, datetime.datetime.now().strftime('%Y%m%d_%H%M%S')))
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w', newline='', encoding='utf-8-sig') as fh:
            w = csv.writer(fh, delimiter=';')
            w.writerow(['тип', 'id', 'состояние', 'вид', 'название', 'напряжение_кВ', 'исходный_слой', 'длина_м'])
            for k, s in sorted(self.last['es'].items()):
                if s in (1, 2):
                    e = self.E[k]
                    w.writerow(['линия', k, STATE_TXT[s], e['kind'], e['name'], e['volt'] or '', e['layer'],
                                round(e['len'], 1)])
            for k, s in sorted(self.last['ns'].items()):
                if s in (1, 2):
                    n = self.N[k]
                    w.writerow(['узел', k, STATE_TXT[s], n['type'], n['name'], '', '', ''])
        return path

    def dead_extent(self):
        fe = [self.E[k]['fid'] for k, s in (self.last or {}).get('es', {}).items() if s in (1, 2)]
        fn = [self.N[k]['fid'] for k in self.off_n]
        r = None
        for lyr, fids in ((self.le, fe), (self.ln, fn)):
            if not fids:
                continue
            rq = QgsFeatureRequest().setFilterFids(fids)
            rq.setSubsetOfAttributes([])
            for f in lyr.getFeatures(rq):
                b = f.geometry().boundingBox()
                if r is None:
                    r = QgsRectangle(b)
                else:
                    r.combineExtentWith(b)
        if r is None:
            return None
        pad = max(r.width(), r.height(), 200) * 0.15
        return QgsRectangle(r.xMinimum() - pad, r.yMinimum() - pad, r.xMaximum() + pad, r.yMaximum() + pad)

    # ---- поиск объекта по клику (через пространственный индекс GeoPackage)
    def pick(self, x, y, tol):
        pt = QgsGeometry.fromPointXY(QgsPointXY(x, y))
        best = None
        rq = QgsFeatureRequest().setFilterRect(QgsRectangle(x - 3 * tol, y - 3 * tol, x + 3 * tol, y + 3 * tol))
        rq.setSubsetOfAttributes([])
        for f in self.ln.getFeatures(rq):
            k = self.fid2nid.get(f.id())
            if k is None:
                continue
            d = f.geometry().distance(pt)
            if self.N[k]['src']:
                d = d / 3.0   # приоритет станций и точек питания (они рисуются крупным значком)
            if d <= tol and (best is None or d < best[0]):
                best = (d, 'node', k)
        if best:
            return best[1], best[2]
        rq = QgsFeatureRequest().setFilterRect(QgsRectangle(x - tol, y - tol, x + tol, y + tol))
        rq.setSubsetOfAttributes([])
        for f in self.le.getFeatures(rq):
            k = self.fid2eid.get(f.id())
            if k is None:
                continue
            d = f.geometry().distance(pt)
            if d <= tol and (best is None or d < best[0]):
                best = (d, 'line', k)
        return (best[1], best[2]) if best else (None, None)

    def find(self, text, limit=500):
        """Поиск линий и станций по наименованию (без учёта регистра). Возвращает (рёбра, узлы)."""
        t = text.strip().lower()
        if not t:
            return [], []
        es = [k for k, e in self.E.items() if t in e['name'].lower()][:limit]
        ns = [k for k, n in self.N.items() if n['name'] and t in n['name'].lower()][:limit]
        return es, ns

    def label(self, kind, k):
        if kind == 'line':
            e = self.E[k]
            v = ('%g кВ' % e['volt']) if e['volt'] else ('шина' if e['kind'] == 'BUS' else e['kind'])
            return 'Линия #%d  %s  (%s, %.0f м)' % (k, e['name'] or '—', v, e['len'])
        n = self.N[k]
        return 'Узел #%d  %s %s' % (k, n['type'], n['name'])

    # ---- производительность
    def benchmark(self, n=None, max_lines=3, seed=1, t0_manual=60.0, t_edge=8.0, t_branch=15.0):
        """Случайные сценарии отключения: время программы vs модель ручного анализа оператором."""
        n = n or (100 if len(self.E) < 20000 else 30)
        rnd = random.Random(seed)
        cand = [k for k, e in self.E.items() if e['kind'] == 'LINE' and self.G.edges[k][0] in self.G.base_nodes]
        rows = []
        for i in range(n):
            off = rnd.sample(cand, rnd.randint(1, max_lines))
            es, ns, dt = self.G.analyze(off, ())
            nd = sum(1 for k, s in es.items() if s == 2 and self.E[k]['kind'] == 'LINE')
            nb = sum(1 for k, s in ns.items() if s == 2 and self.N[k]['type'] == 'BRANCH')
            manual = t0_manual + t_edge * (nd + len(off)) + t_branch * nb
            rows.append((len(off), nd, nb, dt, manual))
        tp = sum(r[3] for r in rows) / n
        tm = sum(r[4] for r in rows) / n
        full = self.run()
        txt = ('Сеть: %d рёбер. Сценариев: %d. Среднее: обесточено %.1f участков, %.1f ответвлений.<br>'
               'Программа: %.1f мс на расчёт (с перерисовкой ~%.0f мс).<br>'
               'Модель ручного анализа: %.0f с (%.1f мин) = %g + %g·N<sub>уч</sub> + %g·N<sub>отв</sub>.<br>'
               '<b>Ускорение ≈ %.0f раз</b> (по полному времени с перерисовкой).'
               % (len(self.E), n, sum(r[1] for r in rows) / n, sum(r[2] for r in rows) / n, tp * 1000,
                  (tp + full['t_total'] - full['t_calc']) * 1000, tm, tm / 60, t0_manual, t_edge, t_branch,
                  tm / max(tp + full['t_total'] - full['t_calc'], 1e-6)))
        return txt, rows

    # ---- ручной замер
    def manual_start(self):
        """Показывает только отключённые оператором элементы (без последствий) и запускает секундомер."""
        es, ns, dt = self.G.analyze(self.off_e, self.off_n)
        self._paint({k: (1 if s == 1 else (3 if s == 3 else 0)) for k, s in es.items()},
                    {k: (1 if s == 1 else (3 if s == 3 else 0)) for k, s in ns.items()})
        self.le.removeSelection()
        iface.setActiveLayer(self.le)
        iface.actionSelect().trigger()
        self.manual_t0 = time.perf_counter()

    def manual_stop(self):
        if self.manual_t0 is None:
            return 'Ручной замер не запущен.'
        tm = time.perf_counter() - self.manual_t0
        self.manual_t0 = None
        sel = {f['edge_id'] for f in self.le.selectedFeatures()}
        r = self.run()
        truth = {k for k, s in r['es'].items() if s == 2 and self.E[k]['kind'] == 'LINE'}
        sel_l = {k for k in sel if self.E[k]['kind'] == 'LINE'}
        tp = len(sel_l & truth)
        fp = len(sel_l - truth)
        fn = len(truth - sel_l)
        self.le.removeSelection()
        rec = tp / max(len(truth), 1)
        txt = ('Ручной метод: <b>%.1f с</b>, выбрано %d участков: верно %d, лишних %d, пропущено %d '
               '(полнота %.0f%%).<br>Программа: <b>%.1f мс</b>, найдено 100%% участков.<br>'
               '<b>Программа быстрее в %.0f раз</b>' % (tm, len(sel_l), tp, fp, fn, 100.0 * rec,
                                                        r['t_total'] * 1000, tm / max(r['t_total'], 1e-6)))
        if 0 < rec < 1:
            txt += (' (при полной ручной обработке ≈ %.0f с, т.е. в %.0f раз)'
                    % (tm / rec, tm / rec / max(r['t_total'], 1e-6)))
        return txt + '.'


# ------------------------------------------------------------------ панель
class OutageDock(QDockWidget):
    def __init__(self, tool):
        super().__init__('Анализ аварийных отключений – %s' % TAG)
        self.setObjectName('OutageDock')
        self.tool = tool
        self.canvas = iface.mapCanvas()
        self.map_tool = PickTool(self.canvas, self.on_pick)
        w = QWidget()
        v = QVBoxLayout(w)
        v.addWidget(QLabel('<b>1. Выбор элементов для отключения</b><br>'
                           '<small>«Выбор на карте»: ЛКМ — отключить/включить линию или узел.<br>'
                           'Или выделите объекты инструментом выделения и нажмите «Добавить выделенные».</small>'))
        h = QHBoxLayout()
        self.ed_find = QLineEdit()
        self.ed_find.setPlaceholderText('Поиск линии / станции по названию…')
        self.ed_find.returnPressed.connect(self.on_find)
        b = QPushButton('Найти')
        b.clicked.connect(self.on_find)
        h.addWidget(self.ed_find)
        h.addWidget(b)
        v.addLayout(h)
        h = QHBoxLayout()
        self.b_pick = QPushButton('Выбор на карте')
        self.b_pick.setCheckable(True)
        self.b_pick.toggled.connect(self.on_pick_mode)
        b_sel = QPushButton('Добавить выделенные')
        b_sel.clicked.connect(self.on_add_selected)
        h.addWidget(self.b_pick)
        h.addWidget(b_sel)
        v.addLayout(h)
        self.lst = QListWidget()
        self.lst.setMaximumHeight(130)
        v.addWidget(self.lst)
        h = QHBoxLayout()
        b = QPushButton('Убрать из списка')
        b.clicked.connect(self.on_remove)
        h.addWidget(b)
        b = QPushButton('Сбросить всё')
        b.clicked.connect(self.on_reset)
        h.addWidget(b)
        v.addLayout(h)
        self.cb_auto = QCheckBox('Пересчитывать автоматически')
        self.cb_auto.setChecked(True)
        v.addWidget(self.cb_auto)
        b = QPushButton('▶  Выполнить анализ отключения')
        f = QFont()
        f.setBold(True)
        b.setFont(f)
        b.setMinimumHeight(30)
        b.clicked.connect(self.on_run)
        v.addWidget(b)
        v.addWidget(QLabel('<b>2. Результат</b>'))
        self.txt = QTextEdit()
        self.txt.setReadOnly(True)
        v.addWidget(self.txt)
        h = QHBoxLayout()
        b = QPushButton('Показать на карте')
        b.clicked.connect(self.on_zoom)
        h.addWidget(b)
        b = QPushButton('Экспорт CSV')
        b.clicked.connect(self.on_export)
        h.addWidget(b)
        v.addLayout(h)
        g = QGroupBox('3. Сравнение с ручным методом')
        gv = QVBoxLayout(g)
        h = QHBoxLayout()
        b = QPushButton('Старт ручного замера')
        b.clicked.connect(self.on_manual_start)
        h.addWidget(b)
        b = QPushButton('Стоп и сравнить')
        b.clicked.connect(self.on_manual_stop)
        h.addWidget(b)
        gv.addLayout(h)
        b = QPushButton('Тест скорости (случайные аварии)')
        b.clicked.connect(self.on_bench)
        gv.addWidget(b)
        v.addWidget(g)
        self.setWidget(w)
        self.refresh_list()

    def refresh_list(self):
        self.lst.clear()
        for k in sorted(self.tool.off_e):
            self.lst.addItem(self.tool.label('line', k))
            self.lst.item(self.lst.count() - 1).setData(Qt.UserRole, ('line', k))
        for k in sorted(self.tool.off_n):
            self.lst.addItem(self.tool.label('node', k))
            self.lst.item(self.lst.count() - 1).setData(Qt.UserRole, ('node', k))

    def after_change(self):
        self.refresh_list()
        if self.cb_auto.isChecked():
            self.on_run()

    def on_pick_mode(self, on):
        if on:
            self.prev_tool = self.canvas.mapTool()
            self.canvas.setMapTool(self.map_tool)
        else:
            self.canvas.unsetMapTool(self.map_tool)
            if getattr(self, 'prev_tool', None):
                self.canvas.setMapTool(self.prev_tool)

    def on_pick(self, pt, tol, right):
        dst = self.tool.le.crs()
        src = self.canvas.mapSettings().destinationCrs()
        if src != dst:
            tr = QgsCoordinateTransform(src, dst, QgsProject.instance())
            p2 = tr.transform(QgsPointXY(pt.x() + tol, pt.y()))
            pt = tr.transform(pt)
            tol = abs(p2.x() - pt.x()) or tol
        if tol > PICK_MAX_TOL:
            iface.messageBar().pushWarning('Анализ отключений', 'Увеличьте масштаб карты для выбора объекта')
            return
        kind, k = self.tool.pick(pt.x(), pt.y(), tol)
        if kind is None:
            iface.messageBar().pushInfo('Анализ отключений', 'Рядом с точкой нет линии или узла сети')
            return
        if kind == 'line':
            self.tool.toggle_line(k)
        else:
            self.tool.toggle_node(k)
        self.after_change()

    def on_find(self):
        es, ns = self.tool.find(self.ed_find.text())
        self.tool.le.selectByIds([self.tool.E[k]['fid'] for k in es])
        self.tool.ln.selectByIds([self.tool.N[k]['fid'] for k in ns])
        if not es and not ns:
            iface.messageBar().pushInfo('Анализ отключений', 'Ничего не найдено')
            return
        box = None
        for lyr in (self.tool.le, self.tool.ln):
            if lyr.selectedFeatureCount():
                b = lyr.boundingBoxOfSelected()
                box = b if box is None else (box.combineExtentWith(b) or box)
        if box is not None:
            box.grow(max(box.width(), box.height(), 500) * 0.1)
            dstc = self.canvas.mapSettings().destinationCrs()
            if self.tool.le.crs() != dstc:
                box = QgsCoordinateTransform(self.tool.le.crs(), dstc, QgsProject.instance()).transformBoundingBox(box)
            self.canvas.setExtent(box)
            self.canvas.refresh()
        self.txt.setHtml('Найдено: линий %d, станций/узлов %d (выделены на карте). '
                         'Нажмите «Добавить выделенные», чтобы отключить их.' % (len(es), len(ns)))

    def on_add_selected(self):
        n = self.tool.add_selected()
        if not n:
            iface.messageBar().pushWarning('Анализ отключений',
                                           'Нет выделенных объектов в слоях «%s» / «%s»' % (TITLE[L_EDGES], TITLE[L_NODES]))
        self.after_change()

    def on_remove(self):
        for it in self.lst.selectedItems():
            kind, k = it.data(Qt.UserRole)
            (self.tool.off_e if kind == 'line' else self.tool.off_n).discard(k)
        self.after_change()

    def on_reset(self):
        self.tool.reset()
        self.refresh_list()
        self.txt.setHtml('Сеть восстановлена в нормальный режим.')

    def on_run(self):
        self.tool.run()
        self.txt.setHtml(self.tool.report())

    def on_zoom(self):
        r = self.tool.dead_extent()
        if r:
            src = self.tool.le.crs()
            dst = self.canvas.mapSettings().destinationCrs()
            if src != dst:
                r = QgsCoordinateTransform(src, dst, QgsProject.instance()).transformBoundingBox(r)
            self.canvas.setExtent(r)
            self.canvas.refresh()

    def on_export(self):
        p = self.tool.export_csv()
        iface.messageBar().pushSuccess('Анализ отключений', 'Отчёт сохранён: ' + p)

    def on_manual_start(self):
        if not self.tool.off_e and not self.tool.off_n:
            QMessageBox.information(self, 'Ручной замер', 'Сначала выберите отключаемые элементы.')
            return
        self.cb_auto.setChecked(False)
        if self.b_pick.isChecked():
            self.b_pick.setChecked(False)
        self.tool.manual_start()
        self.txt.setHtml('<b>Ручной замер запущен.</b><br>Показаны только отключённые элементы (красным). '
                         'Выделите инструментом выделения ВСЕ линии, которые, по-вашему, окажутся обесточены '
                         '(Ctrl — добавить к выделению), затем нажмите «Стоп и сравнить».')

    def on_manual_stop(self):
        self.txt.setHtml(self.tool.manual_stop())
        self.cb_auto.setChecked(True)

    def on_bench(self):
        txt, rows = self.tool.benchmark()
        self.txt.setHtml(txt)

    def closeEvent(self, e):
        if self.b_pick.isChecked():
            self.b_pick.setChecked(False)
        super().closeEvent(e)


def start():
    """Точка входа: строит/загружает сеть и открывает панель."""
    global OUTAGE_TOOL, OUTAGE_DOCK
    prj = QgsProject.instance()
    if REBUILD:
        build_network()
    elif not prj.mapLayersByName(TITLE[L_EDGES]):
        if os.path.exists(OUT_GPKG):
            load_network_layers(OUT_GPKG)
        else:
            build_network()
    else:
        style_layers(*[prj.mapLayersByName(TITLE[n])[0] if prj.mapLayersByName(TITLE[n]) else None
                       for n in (L_EDGES, L_NODES, L_MASTS)])
    OUTAGE_TOOL = OutageTool()
    OUTAGE_TOOL.run()
    try:
        OUTAGE_DOCK.close()
        iface.removeDockWidget(OUTAGE_DOCK)
    except Exception:
        pass
    OUTAGE_DOCK = OutageDock(OUTAGE_TOOL)
    iface.addDockWidget(Qt.RightDockWidgetArea, OUTAGE_DOCK)
    OUTAGE_DOCK.txt.setHtml('Сеть «%s» загружена: %d рёбер, %d узлов, %d источников питания, %d опор.<br>'
                            'Найдите объект по названию или включите «Выбор на карте».'
                            % (TAG, len(OUTAGE_TOOL.E), len(OUTAGE_TOOL.N), len(OUTAGE_TOOL.G.sources), OUTAGE_TOOL.n_masts))
    le = OUTAGE_TOOL.le
    ext = le.extent()
    dst = iface.mapCanvas().mapSettings().destinationCrs()
    if le.crs() != dst:
        ext = QgsCoordinateTransform(le.crs(), dst, prj).transformBoundingBox(ext)
    iface.mapCanvas().setExtent(ext)
    iface.mapCanvas().refresh()
    return OUTAGE_TOOL


OUTAGE_TOOL = globals().get('OUTAGE_TOOL')
OUTAGE_DOCK = globals().get('OUTAGE_DOCK')
start()
print('Готово. Объект OUTAGE_TOOL доступен в консоли: OUTAGE_TOOL.disconnect(lines=[...], nodes=[...]), '
      'OUTAGE_TOOL.find("название"), OUTAGE_TOOL.reset(), OUTAGE_TOOL.benchmark()')
