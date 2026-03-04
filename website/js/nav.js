/* ============================================================
   Geometric Wave Theory — Navigation Component
   Collapsible sections · Mobile hamburger · Active page highlight
   ============================================================ */

(function () {
  'use strict';

  /* --- Navigation Data --- */
  const NAV = [
    { label: 'Home', href: '/index.html', icon: '' },
    { label: 'About', href: '/pages/about.html' },
    { label: 'References', href: '/pages/references.html' },
    {
      label: 'The Lattice',
      items: [
        { label: 'The Medium', href: '/pages/the-lattice.html#medium' },
        { label: 'The Potential V(x)', href: '/pages/the-lattice.html#potential' },
        { label: 'Constant Equivalence', href: '/pages/the-lattice.html#equivalence' },
        { label: "Hooke's Law Completed", href: '/pages/the-lattice.html#hooke' },
      ]
    },
    {
      label: 'Core Derivations',
      items: [
        { label: 'c, \u0127, G from {k, a, \u03b7}', href: '/pages/core-derivations.html#master' },
        { label: '\u03b1 = 1/137.036', href: '/pages/core-derivations.html#alpha' },
        { label: 'Strong coupling \u03b1s', href: '/pages/core-derivations.html#alphas' },
        { label: 'Mass formulas', href: '/pages/core-derivations.html#mass' },
        { label: '\u03a9_\u039b = 2/3', href: '/pages/core-derivations.html#omega' },
      ]
    },
    {
      label: 'Electromagnetism',
      items: [
        { label: 'Overview', href: '/pages/electromagnetism.html' },
        { label: 'Predictions', href: '/pages/electromagnetism.html#predictions' },
      ]
    },
    {
      label: 'Nuclear / QCD',
      items: [
        { label: 'Overview', href: '/pages/nuclear.html' },
        { label: 'Predictions', href: '/pages/nuclear.html#predictions' },
      ]
    },
    {
      label: 'Gravity & Cosmology',
      items: [
        { label: 'Overview', href: '/pages/gravity-cosmology.html' },
        { label: 'Predictions', href: '/pages/gravity-cosmology.html#predictions' },
        { label: 'Nested Well Suppression', href: '/pages/gravity-cosmology.html#nested-wells' },
      ]
    },
    {
      label: 'Particle Physics',
      items: [
        { label: 'Overview', href: '/pages/particle-physics.html' },
        { label: 'Predictions', href: '/pages/particle-physics.html#predictions' },
      ]
    },
    {
      label: 'Quantum Mechanics',
      items: [
        { label: 'Overview', href: '/pages/quantum-mechanics.html' },
        { label: 'Predictions', href: '/pages/quantum-mechanics.html#predictions' },
      ]
    },
    {
      label: 'Multi-Well Cosmology',
      items: [
        { label: 'Parallel Worlds', href: '/pages/multi-well.html#parallel' },
        { label: 'Dark Matter', href: '/pages/multi-well.html#dark-matter' },
      ]
    },
    {
      label: 'Cyclic Cosmology',
      items: [
        { label: 'The Self-Resetting Universe', href: '/pages/cyclic-cosmology.html' },
      ]
    },
    {
      label: 'Resolved Mysteries',
      items: [
        { label: 'Hierarchy Problem', href: '/pages/resolved-mysteries.html#hierarchy' },
        { label: 'Strong CP', href: '/pages/resolved-mysteries.html#strong-cp' },
        { label: 'Information Paradox', href: '/pages/resolved-mysteries.html#information' },
        { label: 'Proton Stability', href: '/pages/resolved-mysteries.html#proton' },
        { label: 'Arrow of Time', href: '/pages/resolved-mysteries.html#arrow' },
        { label: 'Why 3+1 Dimensions', href: '/pages/resolved-mysteries.html#dimensions' },
      ]
    },
    {
      label: 'Engineering',
      items: [
        { label: 'Warp Drive', href: '/pages/engineering.html#warp' },
        { label: 'Lattice Communications', href: '/pages/engineering.html#comms' },
      ]
    },
    {
      label: 'Interactive Tools',
      items: [
        { label: 'Atomic Mass Predictor', href: '/tools/atomic-mass-predictor.html' },
        { label: 'Proton Form Factor', href: '/tools/proton-form-factor.html' },
        { label: 'Wave Atom Modeler', href: '/tools/wave-atom-modeler.html' },
        { label: 'Energy Density Plot', href: '/tools/energy-density-plot.html' },
        { label: 'Cosmology Simulator', href: '/tools/lattice-cosmology-simulator.html' },
        { label: 'Cross-Well Calculator', href: '/tools/cross-well-gravity.html' },
        { label: 'Cyclic Timeline', href: '/tools/cyclic-cosmology.html' },
      ]
    },
    {
      label: 'Calculations',
      items: [
        { label: 'Master Equations', href: '/calculations/calc-master-equations.html' },
        { label: 'Fine Structure \u03b1', href: '/calculations/calc-fine-structure.html' },
        { label: 'Strong Coupling & QCD', href: '/calculations/calc-strong-qcd.html' },
        { label: 'Particle Masses', href: '/calculations/calc-particle-masses.html' },
        { label: 'Electroweak', href: '/calculations/calc-electroweak.html' },
        { label: 'Proton & Nuclear', href: '/calculations/calc-proton-nuclear.html' },
        { label: 'Atomic Physics', href: '/calculations/calc-atomic.html' },
        { label: 'Neutrinos', href: '/calculations/calc-neutrinos.html' },
        { label: 'Cosmology', href: '/calculations/calc-cosmology.html' },
        { label: 'Gravity & GR Tests', href: '/calculations/calc-gravity-gr.html' },
        { label: 'Mixing Angles', href: '/calculations/calc-mixing-angles.html' },
        { label: 'The Lagrangian', href: '/calculations/calc-lagrangian.html' },
      ]
    },
    { label: 'All Predictions', href: '/predictions/index.html' },
  ];

  /* --- Determine base path --- */
  function getBasePath() {
    const path = window.location.pathname;
    // If we're at root or in website/ root
    if (path.endsWith('/index.html') && !path.includes('/pages/') && !path.includes('/tools/') && !path.includes('/predictions/') && !path.includes('/calculations/')) {
      return '.';
    }
    if (path.includes('/pages/') || path.includes('/tools/') || path.includes('/predictions/') || path.includes('/calculations/')) {
      return '..';
    }
    return '.';
  }

  /* --- Resolve href relative to base --- */
  function resolve(href) {
    const base = getBasePath();
    return base + href;
  }

  /* --- Check if current page matches href --- */
  function isActive(href) {
    const current = window.location.pathname;
    const resolved = href.split('#')[0]; // strip hash
    return current.endsWith(resolved) || current.endsWith(resolved.replace(/^\//, ''));
  }

  /* --- Build Navigation DOM --- */
  function buildNav() {
    const sidebar = document.querySelector('.nav-sidebar');
    if (!sidebar) return;

    // Header
    const header = document.createElement('div');
    header.className = 'nav-header';
    header.innerHTML = `<a href="${resolve('/index.html')}">Geometric Wave Theory<span class="nav-subtitle">179 predictions &middot; 0 free parameters</span></a>`;
    sidebar.appendChild(header);

    // Body
    const body = document.createElement('div');
    body.className = 'nav-body';

    NAV.forEach(item => {
      if (item.items) {
        // Section with children
        const section = document.createElement('div');
        section.className = 'nav-section';

        const sectionHeader = document.createElement('div');
        sectionHeader.className = 'nav-section-header';
        sectionHeader.innerHTML = `<span>${item.label}</span><span class="chevron">&#9654;</span>`;
        section.appendChild(sectionHeader);

        const itemsDiv = document.createElement('div');
        itemsDiv.className = 'nav-section-items';

        let sectionHasActive = false;
        item.items.forEach(child => {
          const a = document.createElement('a');
          a.className = 'nav-link';
          a.href = resolve(child.href);
          a.textContent = child.label;
          if (isActive(child.href)) {
            a.classList.add('active');
            sectionHasActive = true;
          }
          itemsDiv.appendChild(a);
        });

        section.appendChild(itemsDiv);

        // Auto-open sections with active links
        if (sectionHasActive) {
          section.classList.add('open');
        }

        // Toggle handler
        sectionHeader.addEventListener('click', () => {
          section.classList.toggle('open');
        });

        body.appendChild(section);
      } else {
        // Top-level link
        const a = document.createElement('a');
        a.className = 'nav-top-link';
        a.href = resolve(item.href);
        a.textContent = item.label;
        if (isActive(item.href)) {
          a.classList.add('active');
        }
        body.appendChild(a);
      }
    });

    sidebar.appendChild(body);
  }

  /* --- Mobile Toggle --- */
  function setupMobile() {
    const toggle = document.querySelector('.nav-toggle');
    const sidebar = document.querySelector('.nav-sidebar');
    const overlay = document.querySelector('.nav-overlay');

    if (!toggle || !sidebar) return;

    toggle.addEventListener('click', () => {
      sidebar.classList.toggle('open');
      if (overlay) overlay.classList.toggle('active');
    });

    if (overlay) {
      overlay.addEventListener('click', () => {
        sidebar.classList.remove('open');
        overlay.classList.remove('active');
      });
    }
  }

  /* --- Init --- */
  document.addEventListener('DOMContentLoaded', () => {
    buildNav();
    setupMobile();
  });
})();
