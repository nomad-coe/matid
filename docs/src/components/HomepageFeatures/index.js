import React from 'react';
import Link from '@docusaurus/Link'
import clsx from 'clsx';
import styles from './styles.module.css';

const FeatureList = [
  {
    title: 'Symmetry-based Clustering',
    link: '/docs/learn/symmetry-based-clustering',
    description: (
      <>
        MatID contains a novel clustering approach called Symmetry-based
        Clustering (SBC) that can be used to identify meaningful components from
        atomistic systems.
      </>
    ),
  },
  {
    title: 'Symmetry Analysis',
    link: '/docs/learn/symmetry-analysis',
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
    linktext: 'Visit NOMAD',
    description: (
      <>
        MatID is being developed at FAIRmat and is powering the automated
        analysis of structures at scale in the NOMAD platform.
      </>
    ),
  },
];

function Feature({title, description, link, linktext}) {
  return (
    <div className={clsx([styles.feature, 'col col--4'])}>
      <div className="text--center padding-horiz--md">
        <h3>{title}</h3>
        <p>{description}</p>
        <Link
          className="button button--secondary button--md"
          to={link}>
            {linktext || 'Learn more'}
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
