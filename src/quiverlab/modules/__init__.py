"""quiverlab.modules: right A-modules over any exact Domain, and their homological
algebra (simples/projectives/injectives, radical/top/socle, Hom/End, minimal
projective resolutions, Ext). Right-module (anti-homomorphism) convention throughout:
m*b = action[b] @ m (m a column), action[x*y] = action[y] @ action[x]. The minimal
resolution engine (resolution.py) is generalized from the bridge obstruction engine
(spec §5 component 7; MIT header there)."""
