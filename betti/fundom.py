from sage.all import QQ

def fundamental_domain(self, group=None):
    """
    Return a fundamental domain for the action of a group of
    symmetries on this polyhedron.

    INPUT:

    - ``group`` -- (optional) if this is given, it must be a finite
      matrix group representing affine or linear transformations on
      the ambient space inducing symmetries of the polyhedron. If
      ``group`` is not given, the
      :meth:`restricted_automorphism_group` is used.
      The transformation matrices act on the left.
      Instead of a group, a list or iterable with the elements is
      also accepted. In that case, it is assumed that those elements
      form a group.
      The group is assumed to act faithfully: no non-trivial
      symmetry should fix every point of the polyhedron.

    OUTPUT: a sub-polyhedron of the given polyhedron which tiles
    the given polyhedron under the action of the group.

    EXAMPLES::

        sage: P = Polyhedron([(0,0), (1,0), (0,1)])
        sage: P.volume()
        1/2
        sage: Q = P.fundamental_domain()
        sage: Q.vertices()
        (A vertex at (0, 0), A vertex at (1/3, 1/3), A vertex at (0, 1/2))
        sage: Q.volume()
        1/12

    A line in 3 dimensions::

        sage: P = Polyhedron([(0,0,0), (2,4,6)])
        sage: P.fundamental_domain().vertices()
        (A vertex at (1, 2, 3), A vertex at (0, 0, 0))

    A non-compact polyhedron::

        sage: P = Polyhedron(lines=[(0,1,0,0)], rays=[(1,0,0,0), (0,0,1,0)])
        sage: Q = P.fundamental_domain()
        sage: Q.Vrepresentation()
        (A line in the direction (0, 1, 0, 0),
         A vertex at (0, 0, 0, 0),
         A ray in the direction (1, 0, 1, 0),
         A ray in the direction (0, 0, 1, 0))

    The 24-cell::

        sage: P = polytopes.twenty_four_cell()
        sage: P.volume()
        2
        sage: Q = P.fundamental_domain(); Q  # long time
        A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 5 vertices
        sage: Q.volume()  # long time
        1/576

    A triangle, first with the full symmetry group, then with a
    subgroup::

        sage: P = Polyhedron([(0,0), (1,0), (0,2)])
        sage: P.volume()
        1
        sage: G = P.restricted_automorphism_group(kind="matrix")
        sage: Q = P.fundamental_domain(G)
        sage: Q.vertices()
        (A vertex at (0, 0), A vertex at (1/3, 2/3), A vertex at (0, 1))
        sage: Q.volume()
        1/6
        sage: H = G.intersection(GL(3,ZZ))
        sage: Q = P.fundamental_domain(H)
        sage: Q.vertices()
        (A vertex at (1, 0), A vertex at (0, 1), A vertex at (0, 0))
        sage: Q.volume()
        1/2

    The 2-dimensional plane with a group which is not a subgroup
    of the restricted automorphism group::

        sage: P = Polyhedron(lines=[(1,0), (0,1)])
        sage: G = MatrixGroup([[[1,-1],[1,0]]])  # C_6
        sage: Q = P.fundamental_domain(G)
        sage: Q.Vrepresentation()
        (A vertex at (0, 0),
         A ray in the direction (1, 2),
         A ray in the direction (-1, 1))

    Now with the trivial group. We get the same polyhedron but
    defined over a field::

        sage: P
        A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 lines
        sage: P.fundamental_domain([identity_matrix(2)])
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines
    """
    # The fundamental domain is rather meaningless for
    # 0-dimensional polyhedra: just return the polyhedron itself.
    if self.dimension() <= 0:
        return self

    if group is None:
        G = self.restricted_automorphism_group(kind="matrixlist")
    else:
        # Convert all group elements to matrices representing affine
        # transformations
        from sage.groups.affine_gps.affine_group import AffineGroup
        AGL = AffineGroup(self.ambient_dim(), self.base_ring().fraction_field())
        G = [AGL(g).matrix() for g in group]

    # We need to find a good "center" point for the fundamental
    # domain. This is a point such that the size of the orbit
    # equals the group order (in other words, the point is not
    # fixed by any non-trivial group element). Ideally, the point
    # is an "easy" point (for example, a vertex or the midpoint of
    # an edge).
    class unique_set(set):
        def add_unique(self, obj):
            "Add ``obj``, which cannot be already in the set"
            if obj in self:
                raise LookupError
            self.add(obj)

    seen = unique_set()
    def is_good_center_point(pt):
        try:
            for g in G:
                gpt = g * pt
                gpt.set_immutable()
                seen.add_unique(gpt)
            return True
        except LookupError:
            return False

    # First, try all V-representation objects.
    V = [v.homogeneous_vector() for v in self.Vrepresentation()]
    for center in V:
        if is_good_center_point(center):
            break
    else:
        # Try the midpoint of all edges.
        for f in self.faces(1):
            v = f.ambient_Vrepresentation(0).homogeneous_vector()
            w = f.ambient_Vrepresentation(1).homogeneous_vector()
            # homogeneous coordinates, no need to divide by 2
            center = v + w
            if is_good_center_point(center):
                break
        else:
            # Finally, try random points by summing randomly chosen
            # V-representation objects.
            from sage.misc.randstate import random
            l = len(V)

            # Start with some point which we not have already
            # tried.
            center = V[0] + V[0] + V[1]
            while not is_good_center_point(center):
                center += V[random() % l]

    # This matrix defines an inner product <v,w> = v^t M w chosen
    # such that <v,w> = <g(v),g(w)> for all g in G.
    M = sum(g.transpose() * g for g in G)

    # A fundamental domain is given by the points of the given
    # polyhedron which are closer to "center" than any other point
    # in the orbit of "center". To define "closer", we use the
    # inner product defined by M. We start from the equations
    # and inequalities defining the polyhedron itself.
    eqns = list(self.equations())
    ieqs = list(self.inequalities())
    for g in G:
        w = g * center
        if w == center:
            continue
        # Let m be the middle point (c+w)/2. We need to find all x
        # such that <c-w, x-m> >= 0. We rewrite this as
        # <2(c-w), x> - <c-w, c+w> >= 0.
        diff = center - w
        bra = diff * M  # bra * x = <c-w, x>

        # Verify that <c-w, c-w> != 0, otherwise the inequality is
        # trivial.
        assert bra * diff != 0

        left = bra + bra
        const = left[-1] - bra * (center + w)
        ieqs.append([const] + list(left)[:-1])

    parent = self.parent().base_extend(QQ)
    return parent.element_class(parent, None, [ieqs, eqns])
