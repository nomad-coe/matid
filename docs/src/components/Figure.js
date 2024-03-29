import React from 'react'
import useBaseUrl from '@docusaurus/useBaseUrl'

export default function Figure({ src, caption }) {
  return (
    <figure style={{padding: '0.5rem', maxWidth: "36rem", marginRight: "auto", marginLeft: "auto"}}>
      <img src={useBaseUrl(src)} alt={caption} />
      <figcaption>{`Figure: ${caption}`}</figcaption>
    </figure>
  )
}