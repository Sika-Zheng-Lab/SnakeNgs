"""Add dependency-level pulses to a snakevision SVG without changing its graph."""
from graphlib import TopologicalSorter
from pathlib import Path
import re
import xml.etree.ElementTree as ET

SVG = 'http://www.w3.org/2000/svg'
ET.register_namespace('', SVG)


def animate(source, destination, dot):
    labels = dict(re.findall(r'(\d+)\[label\s*=\s*"([^"]+)"', dot))
    parents = {node: set() for node in labels}
    for parent, child in re.findall(r'(\d+)\s*->\s*(\d+)', dot):
        parents[child].add(parent)
    ranks = {}
    for node in TopologicalSorter(parents).static_order():
        ranks[node] = max((ranks[p] + 1 for p in parents[node]), default=0)
    by_name = {labels[node]: rank for node, rank in ranks.items()}
    tree = ET.parse(source)
    root = tree.getroot()
    circles = [node for node in root.iter(f'{{{SVG}}}circle') if node.get('id', '').startswith('N')]
    if not circles:
        raise ValueError(f'No rule nodes in {source}')
    duration = (max(by_name[node.get('id')[1:]] for node in circles) + 2) * 0.8
    peak, end = 100 * 0.35 / duration, 100 * 0.7 / duration
    root.set('tabindex', '0')
    root.set('role', 'img')
    root.set('aria-labelledby', 'flow-title flow-description')
    title = ET.Element(f'{{{SVG}}}title', id='flow-title')
    title.text = f'{Path(source).stem}: animated rulegraph'
    description = ET.Element(f'{{{SVG}}}desc', id='flow-description')
    description.text = ('Nodes pulse in dependency order; parallel branches may pulse together. '
                        'This illustrates dependencies, not execution time or live progress. '
                        'Hover a node or focus this SVG to pause. Reduced-motion preferences disable animation.')
    style = ET.Element(f'{{{SVG}}}style')
    style.text = f'''
.flow-node {{ transform-box: fill-box; transform-origin: center; animation: rule-pulse {duration:.2f}s ease-in-out infinite; }}
@keyframes rule-pulse {{
  0%, {end:.4f}%, 100% {{ transform: scale(1); }}
  {peak:.4f}% {{ transform: scale(1.45); }}
}}
svg:has(.flow-node:hover) .flow-node, svg:focus .flow-node {{ animation-play-state: paused; }}
@media (prefers-reduced-motion: reduce) {{ .flow-node {{ animation: none; }} }}
'''
    for circle in circles:
        rank = by_name[circle.get('id')[1:]]
        circle.set('class', 'flow-node')
        circle.set('style', f'animation-delay: {rank * 0.8:.2f}s')
        circle.set('data-level', str(rank))
    for index, node in enumerate([title, description, style]):
        root.insert(index, node)
    tree.write(destination, encoding='unicode')
