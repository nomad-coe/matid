import React from 'react';
import Link from '@docusaurus/Link'
import clsx from 'clsx';
import styles from './styles.module.css';

const FeatureList = [
  {
    title: 'Identify components using translational symmetry',
    link: '/docs/clustering/clustering-basics',
    description: (
      <>
        There are several tools for building atomistic systems from components,
        but the reverse process of identifying these components from an existing
        system is much harder. This is one of the main tasks for
        which MatID was designed for.
      </>
    ),
  },
  {
    title: 'Symmetry Analysis',
    link: '/docs/symmetry/symmetry-basics',
    description: (
      <>
        MatID contains several symmetry routines for analyzing structures. The
        basic features are based on the excellent {<a
        href="https://github.com/spglib/spglib">spglib library</a>}, but with an
        extended feature set that makes MatID unique.
      </>
    ),
  },
  {
    title: 'Powering identification in NOMAD',
    link: 'https://nomad-lab.eu/nomad-lab/',
    description: (
      <>
        MatID is being developed at FAIRmat and is powering the automated
        analysis of structures at scale in the NOMAD platform.
      </>
    ),
  },
];

function Feature({Svg, title, description, link}) {
  return (
    <div className={clsx([styles.feature, 'col col--4'])}>
      <div className="text--center padding-horiz--md">
        <h3>{title}</h3>
        <p>{description}</p>
        <Link
          className="button button--secondary button--md"
          to={link}>
          Learn more
        </Link>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
