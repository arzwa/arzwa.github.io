---
layout: post
title: Basics of category theory
---

There are five basic concepts in category theory:

- Category
- Functor
- Natural transformations
- Universality
- Adjoints

## Categories

A **category** $C$ is defined by

- A class $Obj(C)$ of objects where $A \in C$ means $A \in Obj(C)$. (i)
- For each pair $A,B \in C$ $\exists$ a set $hom_C(A,B)$ for $(A,B)$. (ii)
- Distinct hom-sets are disjoint (not *really* necessary). (iii)

The elements of $hom_C(A,B)$ called *morphisms, maps* or *arrows*. The class of all morphisms of $C$ is referred to as $Mor(C)$. These are not to be confused with functions in the ordinary sense. Despite this, often functional notation is used to specify morphisms, like $f: A \rightarrow B$. $A$ in this example is called the domain of morphism $f$, while $B$ is called the co-domain. Again do not confuse the meaning of these terms with their counterparts in the case of functions. In general $hom(A,B) \neq hom(B,A)$. A category is also defined by the property of *compositionality*:

- If $f: A \rightarrow B$, $g: B \rightarrow C$, then $\exists$ $g \circ f : A \rightarrow C$. This composition has to be associative, meaning that $f \circ (g \circ h) = (f \circ g) \circ h$. (iv)

Finally a category also has *identity morphisms*:

- For each $A \in C$ $\exists id_A: A \rightarrow A$ such that for $f: A \rightarrow B$, $id_B \circ f = f$ and $f \circ id_A = f$. (v)

If $Obj(C)$ and $Mor(C)$ are sets, than the category is said to be *small*. If the hom-classes are actually sets, then the category is *locally small* (above definition of a category is thus the definition of locally small categories).

As an example of a category consider the category of **sets**, which has as objects all sets, as morphisms set functions and as hom-set $hom(A,B)$ the set of all set functions from $A$ to $B$. A similar example is the category of all relations $R$, where the objects are all sets and $hom(R,S)$ is the set of all relations from $R$ to $S$.

## Functors

**Functors** are morphisms between categories. To establish such a morphism between two categories, $C$ and $D$, we need two maps, namely

$$ F: Obj(C) \rightarrow Obj(D) $$

and

$$ F: Mor(C) \rightarrow Mor(D) $$

where in the *covariant* case (a covariant functor)

$$ F: hom_C(A,B) \rightarrow hom_D(FA, FB) $$

and in the *contravariant* case

$$ F: hom_C(A,B) \rightarrow hom_D(FB, FA) $$
