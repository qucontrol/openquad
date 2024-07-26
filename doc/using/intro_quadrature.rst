.. _quadrature-basics:

What is quadrature?
-------------------

The basic idea of quadrature is to obtain the numerical approximation for the
integral of the function :math:`f(x)` over the *domain* :math:`\mathcal{D}` in
terms of a weighted sum over discrete points,

.. math::

   \int_{\mathcal{D}} f(x) \mathrm{d}x \approx \sum_{i=1}^n f(x_i) w_i
   \;.

Different approaches exist to determine the optimal number and distribution of
*sample points*, :math:`x_i`, and their corresponding *weights*, :math:`w_i`,
for a given domain. This gives rise to a variety of classes of quadrature
methods. See our :ref:`advanced guide <background>` for a brief
overview of different methods.

Based on `Fubini's theorem`_, lower-dimensional quadratures can be combined to
construct tensor product rules for muti-dimensional integrals.

.. _Fubini's theorem: https://en.wikipedia.org/wiki/Fubini%27s_theorem
