import React, { useState, useCallback } from 'react'
import clsx from 'clsx'
import Link from '@docusaurus/Link'
import useDocusaurusContext from '@docusaurus/useDocusaurusContext'
import Layout from '@theme/Layout'
import HomepageFeatures from '@site/src/components/HomepageFeatures'

import styles from './index.module.css'
import Logo from '@site/static/img/logo.png'

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
 * Displays logo, motivational message and get started link.
 */
function GetStarted() {
  const {siteConfig} = useDocusaurusContext()

  return <div className={styles.column}>
    <img src={Logo} className={styles.logo}/>
    <p className="hero__subtitle">{siteConfig.tagline}</p>
    <div className={styles.buttons}>
      <Link
        className="button button--secondary button--lg"
        to="/docs/introduction">
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
    className={styles.stackContainer}
  >
    <img
      onMouseMove={handleMouseMove}
      onMouseLeave={handleMouseLeave}
      src={require(`@site/static/img/${stack}_cropped.png`).default}
      className={styles.stack}
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
