import React, { useState, useCallback } from 'react'
import clsx from 'clsx'
import Link from '@docusaurus/Link'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import Layout from '@theme/Layout'
import HomepageFeatures from '@site/src/components/HomepageFeatures'

import styles from './index.module.css'
import logo from '@site/static/img/logo.png'
import all from '@site/static/img/all_cropped.png'
import cu from '@site/static/img/cu_cropped.png'
import graphene from '@site/static/img/graphene_cropped.png'
import mos2 from '@site/static/img/mos2_cropped.png'

function HomepageHeader() {
  const {siteConfig} = useDocusaurusContext()
  return (
    <header className={clsx('hero hero--primary', styles.heroBanner)}>
      <div className="container">
        <div className={styles.row}>
          <div className={clsx(styles.leftCol)}>
            <GetStarted />
          </div>
          <div className={clsx(styles.rightCol)}>
            <Stack />
          </div>
        </div>
      </div>
    </header>
  );
}

/**
 * Displays logo, motivational message and ed link.
 */
function GetStarted() {
  const {siteConfig} = useDocusaurusContext()

  return <div className={styles.column}>
    <img src={logo.src} className={styles.logo} />
    <p className="hero__subtitle">{siteConfig.tagline}</p>
    <div className={styles.buttons}>
      <Link
        className="button button--secondary button--lg"
        to="/docs/get-started/quick-start">
        Get started
      </Link>
    </div>
  </div>
}

/**
 * Displays an interative image of a clustered system.
 */
function Stack() {
  const [stack, setStack] = useState('all')

  const handleMouseMove = useCallback((event) => {
    const imageHeight = event.target.height
    const y = event.nativeEvent.offsetY
    const fraction = y / imageHeight
    let stack
    if (fraction < 0.25) {
      stack = 'graphene'
    } else if (fraction < 0.57) {
      stack = 'mos2'
    } else {
      stack = 'cu'
    }
    setStack(stack)
  }, [])

  const handleMouseLeave = useCallback(() => setStack('all'), [])

  return <div
    onMouseMove={handleMouseMove}
    onMouseLeave={handleMouseLeave}
    className={styles.stackContainer}
  >
    <img
      src={all.src}
      className={styles.stack}
      style={{visibility: stack === 'all' ? 'visible' : 'hidden', position: 'absolute'}}
    />
    <img
      src={cu.src}
      className={styles.stack}
      style={{visibility: stack === 'cu' ? 'visible' : 'hidden', position: 'absolute'}}
    />
    <img
      src={mos2.src}
      className={styles.stack}
      style={{visibility: stack === 'mos2' ? 'visible' : 'hidden', position: 'absolute'}}
    />
    <img
      src={graphene.src}
      className={styles.stack}
      style={{visibility: stack === 'graphene' ? 'visible' : 'hidden', position: 'absolute'}}
    />
  </div>
}

export default function Home() {
  const {siteConfig} = useDocusaurusContext()
  return (
    <Layout
      title={`${siteConfig.title}`}
      description={`${siteConfig.tagline}`}>
      <HomepageHeader />
      <main>
        <HomepageFeatures />
      </main>
    </Layout>
  );
}
