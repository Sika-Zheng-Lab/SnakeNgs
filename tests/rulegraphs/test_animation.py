import sys
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from animate_rulegraph import animate, SVG


def test_parallel_dependency_order_and_accessibility(tmp_path):
    source = tmp_path / 'source.svg'
    source.write_text(f'<svg xmlns="{SVG}"><circle id="Nfirst"/><circle id="Nleft"/>'
                      '<circle id="Nright"/><circle id="Njoin"/><text>first</text></svg>')
    dot = ('0[label = "first"]; 1[label = "left"]; 2[label = "right"]; '
           '3[label = "join"]; 0 -> 1; 0 -> 2; 1 -> 3; 2 -> 3;')
    output = tmp_path / 'animated.svg'
    animate(source, output, dot)
    tree = ET.parse(output).getroot()
    levels = {n.get('id'): n.get('data-level') for n in tree.iter(f'{{{SVG}}}circle')}
    assert levels == {'Nfirst': '0', 'Nleft': '1', 'Nright': '1', 'Njoin': '2'}
    assert tree.find(f'{{{SVG}}}text').text == 'first'
    css = tree.find(f'{{{SVG}}}style').text
    assert 'prefers-reduced-motion: reduce' in css
    assert 'animation: none' in css
    assert 'animation-play-state: paused' in css
    assert tree.get('tabindex') == '0'
    # Source stays static and can be used as the reduced-motion alternative.
    assert 'animation' not in source.read_text()
